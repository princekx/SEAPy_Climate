"""Compare saved surge and wave products and build the focused assessment."""

import json

import matplotlib.pyplot as plt
import numpy as np
from cf_units import Unit

import analyse_outputs as base


ROLES = tuple(base.RUNS)
COLORS = {"control": "#2471a3", "experiment": "#c0392b", "obs": "#222222"}
LABELS = {"control": "GC5 control", "experiment": "GC6 experiment", "obs": "Observations"}
LAT = np.arange(-10, 25.1, 2.5)
LON = np.arange(95, 130.1, 2.5)
WAVE_LAT = np.arange(-10, 10.1, 2.5)
WAVE_LON = np.arange(100, 130.1, 2.5)


def load(group, role, suffix):
    runid = base.RUNS[role][0]
    return base.load_field(base.folder(group, role) / f"{runid}_{suffix}.nc")


def spatial(group, role, suffix, latitudes=LAT, longitudes=LON):
    return base.sample(load(group, role, suffix), latitudes, longitudes)


def fractions(role, variable, wave, season):
    filtered = load("SEAPy_output", role, f"{variable}_{wave}_{season}_variance")
    total = load("SEAPy_output", role, f"{variable}_ALL_{season}_variance")
    assert filtered["dims"] == total["dims"]
    assert Unit(filtered["units"]) == Unit(total["units"])
    assert all(np.array_equal(filtered["coords"][dim], total["coords"][dim])
               for dim in filtered["dims"])
    threshold = 0.01 * (np.nanmax(total["data"]) - np.nanmin(total["data"]))
    ratio = dict(filtered)
    ratio["data"] = np.divide(100 * filtered["data"], total["data"],
                              out=np.full_like(filtered["data"], np.nan),
                              where=total["data"] > max(threshold, 0))
    return base.sample(ratio, WAVE_LAT, WAVE_LON)


def hovmoller(role, variable, wave):
    field = load("SEAPy_output", role, f"{variable}_{wave}_refEquator_composite")
    assert field["data"].shape[0] == 101
    latitudes = np.arange(-5 if wave == "Kelvin" else -10,
                          5.1 if wave == "Kelvin" else 10.1, 2.5)
    longitudes = np.arange(40, 180.1, 2.5)
    sampled = base.sample(field, latitudes, longitudes)
    values = np.nanmean(sampled, axis=-2)
    assert np.isfinite(values).all() and np.nanmax(np.abs(values)) < 1000
    return np.arange(-50, 51), longitudes, values


def spectrum(role, variable):
    field = load("MJOPy_output", role, f"{variable}_Ratio_Spec_Sym")
    assert field["dims"] == ("frequency", "wavenumber")
    frequencies = field["coords"]["frequency"]
    waves = field["coords"]["wavenumber"]
    selected_frequency = (frequencies >= 1 / 96) & (frequencies <= 0.3)
    selected_wave = (waves >= -15) & (waves <= 15)
    return frequencies[selected_frequency], waves[selected_wave], field["data"][np.ix_(selected_frequency, selected_wave)]


def evidence():
    result = {"surge_anomalies": {}, "fractions": {}, "wave_composites": {}, "spectra": {}}
    weights = np.cos(np.deg2rad(LAT))[:, None]
    for event in ("CS", "CES", "ES", "MS"):
        result["surge_anomalies"][event] = {}
        for variable in ("precipitation_flux", "x_wind_850hPa", "y_wind_850hPa"):
            arrays = [spatial("SEAPy_output", role, f"{event}_{variable}_composite") -
                      spatial("SEAPy_output", role, f"NDJF_{variable}_composite") for role in ROLES]
            result["surge_anomalies"][event][variable] = base.compare_arrays(arrays, weights)
    for variable in ("precipitation_flux", "toa_outgoing_longwave_flux"):
        result["fractions"][variable] = {}
        result["wave_composites"][variable] = {}
        for wave in ("MJO", "Kelvin", "Rossby"):
            for season in ("ndjf", "jjas"):
                arrays = [fractions(role, variable, wave, season) for role in ROLES]
                result["fractions"][variable][f"{wave}_{season}"] = base.compare_arrays(
                    arrays, np.cos(np.deg2rad(WAVE_LAT))[:, None])
            arrays = [hovmoller(role, variable, wave)[2] for role in ROLES]
            lag_limit = 10 if wave == "Kelvin" else 35
            selected = np.abs(np.arange(-50, 51)) <= lag_limit
            result["wave_composites"][variable][wave] = base.compare_arrays(
                [array[selected] for array in arrays], np.ones((1, arrays[0].shape[1])))
        spectra = [spectrum(role, variable) for role in ROLES]
        assert all(np.array_equal(item[0], spectra[0][0]) and np.array_equal(item[1], spectra[0][1])
                   for item in spectra)
        result["spectra"][variable] = base.compare_arrays(
            [item[2] for item in spectra], np.ones_like(spectra[0][2]))
    (base.HERE / "focused_evidence.json").write_text(json.dumps(result, indent=2))
    return result


def surge_figure(result):
    from PIL import Image, ImageChops, ImageDraw, ImageFont

    width, height = 1100, 800
    canvas = Image.new("RGB", (3 * width, 2 * height), "white")
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype("DejaVuSans.ttf", 23)
    for row_index, event in enumerate(("CS", "CES")):
        for column, role in enumerate(ROLES):
            runid = base.RUNS[role][0]
            path = base.folder("SEAPy_output", role) / f"{runid}_{event}_precip_winds850_composite.png"
            with Image.open(path) as source:
                image = source.convert("RGB")
                difference = ImageChops.difference(image, Image.new("RGB", image.size, "white"))
                bounds = difference.convert("L").point(lambda value: 255 if value > 12 else 0).getbbox()
                if bounds is not None:
                    image = image.crop((max(bounds[0] - 8, 0), max(bounds[1] - 8, 0),
                                        min(bounds[2] + 8, image.width), min(bounds[3] + 8, image.height)))
                image.thumbnail((width - 16, height - 45), Image.Resampling.LANCZOS)
                canvas.paste(image, (column * width + (width - image.width) // 2,
                                    row_index * height + 40 + (height - 45 - image.height) // 2))
            draw.text((column * width + 15, row_index * height + 8),
                      f"{'abcdef'[row_index * 3 + column]}) {event}: {LABELS[role]}",
                      fill="black", font=font)
    canvas.save(base.HERE / "focused_figure1_surges.png")


def overlay(ax, longitude, lag, arrays, level):
    artist = ax.contourf(longitude, lag, arrays[2], levels=np.linspace(-2 * level, 2 * level, 17),
                         cmap="RdBu_r", extend="both")
    for role, values in zip(ROLES, arrays):
        styles = {"control": "solid", "experiment": "dashed", "obs": "dotted"}
        contours = ax.contour(longitude, lag, values, levels=[-level, level], colors=[COLORS[role]],
                              linewidths=2.2, linestyles=styles[role])
        ax.clabel(contours, fmt=lambda value: f"{value:+g}", fontsize=8, inline=True)
    ax.set_xlabel("Longitude (E)")
    ax.set_ylabel("Lag (days)")
    return artist


def variance_bars(ax, result, variable):
    keys = [f"{wave}_{season}" for season in ("ndjf", "jjas") for wave in ("MJO", "Kelvin", "Rossby")]
    positions = np.arange(len(keys))
    for offset, role in zip((-0.25, 0, 0.25), ROLES):
        values = [result["fractions"][variable][key][f"{role}_mean"] for key in keys]
        ax.bar(positions + offset, values, width=0.24, color=COLORS[role], label=LABELS[role])
    ax.set_xticks(positions, ["MJO\nNDJF", "Kel\nNDJF", "ER\nNDJF", "MJO\nJJAS", "Kel\nJJAS", "ER\nJJAS"])
    ax.set_ylabel("Mean local variance fraction (%)")
    ax.grid(axis="y", alpha=0.2)


def wave_figure(result):
    fig, axes = plt.subplots(3, 2, figsize=(10.5, 7.5), constrained_layout=True)
    for column, variable in enumerate(("precipitation_flux", "toa_outgoing_longwave_flux")):
        ax = axes[0, column]
        for role in ROLES:
            frequency, waves, values = spectrum(role, variable)
            if role == "obs":
                ax.contourf(waves, frequency, values, levels=np.arange(0.5, 2.76, 0.25),
                            cmap="Greys", extend="both", alpha=0.6)
            ax.contour(waves, frequency, values, levels=[1.2, 1.5], colors=[COLORS[role]],
                       linewidths=1.6, linestyles=["dotted", "solid"], zorder=3)
        ax.set_title(f"{'ab'[column]}) {'Rainfall' if column == 0 else 'OLR'} symmetric spectrum")
        ax.set_xlabel("Zonal wavenumber (positive: eastward)")
        ax.set_ylabel("Frequency (cycles/day)")
    for column, wave in enumerate(("Kelvin", "MJO")):
        lag, longitude, _ = hovmoller("obs", "toa_outgoing_longwave_flux", wave)
        arrays = [hovmoller(role, "toa_outgoing_longwave_flux", wave)[2] for role in ROLES]
        limit = 10 if wave == "Kelvin" else 35
        selected = np.abs(lag) <= limit
        level = 2 if wave == "Kelvin" else 3
        artist = overlay(axes[1, column], longitude, lag[selected],
                         [array[selected] for array in arrays], level)
        axes[1, column].set_title(f"{'cd'[column]}) {wave}: OLR propagation")
        fig.colorbar(artist, ax=axes[1, column], label="Observed OLR anomaly (W/m2)", pad=0.02)
    variance_bars(axes[2, 0], result, "precipitation_flux")
    axes[2, 0].set_title("e) Rainfall variance contribution")
    variance_bars(axes[2, 1], result, "toa_outgoing_longwave_flux")
    axes[2, 1].set_title("f) OLR variance contribution")
    from matplotlib.lines import Line2D
    styles = {"control": "solid", "experiment": "dashed", "obs": "dotted"}
    fig.legend(handles=[Line2D([0], [0], color=COLORS[role], linewidth=2.5,
                               linestyle=styles[role], label=LABELS[role]) for role in ROLES],
               loc="outside lower center", ncol=3, frameon=False)
    fig.savefig(base.HERE / "focused_figure2_waves.png", dpi=240)
    plt.close(fig)


def self_test():
    base.self_test()
    test = base.compare_arrays([np.array([[2.]]), np.array([[1.]]), np.array([[0.]])], np.ones((1, 1)))
    assert test["rmse_change_percent"] == -50
    assert Unit("W/m^2").is_convertible(Unit("W m-2"))
    print("PASS: focused comparison and unit checks")


def document(result):
    from docx import Document
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Inches, Pt

    report = Document()
    section = report.sections[0]
    section.page_width = Inches(8.27)
    section.page_height = Inches(11.69)
    section.top_margin = Inches(0.55)
    section.bottom_margin = Inches(0.55)
    section.left_margin = Inches(0.6)
    section.right_margin = Inches(0.6)
    style = report.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(10.5)
    style.paragraph_format.space_after = Pt(5)
    style.paragraph_format.line_spacing = 1.0
    for heading in ("Title", "Heading 1", "Heading 2"):
        report.styles[heading].font.name = "Calibri"
    report.styles["Title"].font.size = Pt(18)
    report.styles["Heading 1"].font.size = Pt(13)
    report.styles["Heading 2"].font.size = Pt(11)
    footer = section.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    footer.add_run("GC6 MJO and SEA assessment | 30 October 2026 | ")
    page_number = OxmlElement("w:fldSimple")
    page_number.set(qn("w:instr"), "PAGE")
    footer._p.append(page_number)

    def field(label, text):
        paragraph = report.add_paragraph()
        paragraph.add_run(label + " ").bold = True
        paragraph.add_run(text)

    def paragraph(title, text):
        field(title, text)

    def caption(text):
        item = report.add_paragraph(text)
        item.paragraph_format.space_after = Pt(4)
        for run in item.runs:
            run.font.size = Pt(8)

    def figure(name):
        item = report.add_paragraph()
        item.paragraph_format.space_after = Pt(2)
        item.add_run().add_picture(str(base.HERE / name), width=Inches(7.0))

    report.add_heading("Template for Assessment of GC6 Standard Tests", 0)
    report.add_paragraph("Focused assessment of existing MJO and Southeast Asian diagnostic outputs")
    field("Assessment date / requested response:", "30 October 2026 / 3 November 2026.")
    field("Assessment area:", "Cold-surge rainfall and circulation; equatorial-wave propagation; fractional precipitation and OLR variance; wavenumber-frequency spectra; MJO-filtered lag correlations. Keywords: northeast monsoon, Maritime Continent, intraseasonal variability.")
    field("Contributor(s) and organisation(s):", "To be supplied. Dataset-provider acknowledgement statements require confirmation before publication.")
    field("Experiment(s):", "GC6.0 LFRic u-ea703, GC6.0 CPL C144; saved model record December 1981-November 2008.")
    field("Type of experiments:", "Coupled climate simulations, according to run labels. Forcing protocol, ensemble membership and spin-up treatment are not established by the saved outputs.")
    field("Control simulations:", "GC5.0 UM u-cv318, GC5.0 CPL C144; same saved model period.")
    field("Observational datasets:", "Existing obs_ERA5_GPM reference, labelled ERA5/GPM, with saved RMM/event records spanning 2006-2022. Variable-specific product versions and OLR provenance remain unconfirmed. No additional reference data introduced.")
    field("Tools / preprocessing:", "Read-only analysis of MJOPy and SEAPy NetCDF fields and PNGs. Comparisons use a common 2.5-degree grid and common finite masks. Surge anomalies subtract each dataset's NDJF climatology; SEA scores cover 95-130E, 10S-25N. Wave-fraction scores cover 100-130E, 10S-10N, with cosine-latitude weights. Fractions are area means of local filtered/total variance ratios, not ratios of regional means. Spectra use native shared coordinates. See captions for propagation support.")
    field("Location of data used for plots:", "/data/users/prince.xavier/SEAPy_climate/{MJOPy_output,SEAPy_output}/, within u-cv318_GC5.0_CPL_C144/, u-ea703_GC6.0_CPL_C144/ and obs_ERA5_GPM/. Full paths: source_inventory.csv. A verified centrally archived MASS URI has not been supplied.")
    field("One-line summary:", "GC6 improves most wave-variance fractions and selected OLR measures, but worsens cold-surge rainfall/wind anomalies and shows no consistent improvement in rainfall propagation or eastward spectral organisation.")
    field("Any additional information/comments:", "Results are descriptive, not significance-tested: periods differ and internal variability is unquantified. CS dates outside NDJF number 132/1,052 (control), 213/1,255 (experiment) and 76/693 (reference); event-minus-NDJF anomalies are therefore seasonally contaminated. Rainfall-variance metadata are incompatible across datasets: absolute variance comparisons are withheld, while matching units within each filtered/total pair cancel in fractions. The source wave-composite routine allocates an uninitialised array and skips edge peaks without removing their slots; finite, plausible saved values do not exclude contamination. Composite changes are provisional pending an edge-safe regeneration. The package's MJO correlation plot titles say OLR reference, but source code uses each variable's own 80-100E, 10S-10N reference; filtering provenance requires confirmation. No propagation speed or causal attribution is claimed.")
    field("Supporting evidence:", "focused_evidence.json; scalar_metrics.csv; sampling_and_events.json. Reproduce with build_focused_report.py. References: Wheeler & Kiladis (1999), JAS 56, 374-399; Wheeler & Hendon (2004), MWR 132, 1917-1932.")

    report.add_page_break()
    report.add_heading("Detailed Summary (1 of 2): Surge Impacts and Variance", 1)
    figure("focused_figure1_surges.png")
    caption("Figure 1. Existing package plots of event-minus-own-NDJF rainfall (shading, mm/day) and 850-hPa wind (vectors, m/s), reproduced without redrawing their data. CS above, CES below; columns: GC5, GC6, reference. Original grids, colour scales, wind keys and titles are retained. The saved event samples are not restricted to NDJF. Numerical errors in the text use the separate common-grid anomaly comparison.")
    surge = result["surge_anomalies"]
    cs = surge["CS"]["precipitation_flux"]
    ces = surge["CES"]["precipitation_flux"]
    paragraph("Cold-surge response deteriorates.",
              f"CS rainfall-anomaly RMSE rises from {cs['control_rmse']:.2f} to {cs['experiment_rmse']:.2f} mm/day (+{cs['rmse_change_percent']:.0f}%), and pattern correlation falls from {cs['control_correlation']:.2f} to {cs['experiment_correlation']:.2f}. Both models miss the observed wet maximum over the southern South China Sea/northern Borneo and place too much enhancement south of the equator. CS meridional-wind anomaly RMSE increases 18%. Thus the earlier improvement in absolute CS rainfall is partly a background-state gain, not an improved event response.")
    paragraph("Cross-equatorial surges are mixed.",
              f"CES rainfall-anomaly RMSE improves {abs(ces['rmse_change_percent']):.0f}% ({ces['control_rmse']:.2f} to {ces['experiment_rmse']:.2f} mm/day), but the U850/V850 anomaly errors increase 17%/19%. The wet southern Maritime Continent response remains excessive and the dry response is spatially displaced. ES and MS rainfall-anomaly errors also rise, by 21% and 12%; better event frequency or duration does not establish a better circulation-rainfall response.")
    rain = result["fractions"]["precipitation_flux"]
    winter = rain["MJO_ndjf"]
    summer = rain["MJO_jjas"]
    paragraph("Fractional rainfall variance improves selectively (Figure 2e).",
              f"NDJF Kelvin contribution increases from 5.78% to 6.44% (reference 6.74%), and JJAS Rossby excess decreases from 6.31% to 5.65% (3.87%). Local-fraction spatial RMSE improves for Kelvin (14% NDJF; 8% JJAS) and Rossby (26%; 22%). However, the MJO contribution falls from {winter['control_mean']:.2f}% to {winter['experiment_mean']:.2f}% in NDJF (reference {winter['obs_mean']:.2f}%) and from {summer['control_mean']:.2f}% to {summer['experiment_mean']:.2f}% in JJAS ({summer['obs_mean']:.2f}%). NDJF MJO spatial RMSE improves 16% through redistribution despite a worse mean deficit; JJAS RMSE worsens 8%. Fractional gains are not proof of improved absolute rainfall variance, and wave fractions are not additive because filter bands can overlap.")

    report.add_page_break()
    report.add_heading("Detailed Summary (2 of 2): Spectra and Propagation", 1)
    figure("focused_figure2_waves.png")
    caption("Figure 2. Blue/red/black: GC5/GC6/reference. (a,b) Full-record symmetric power/background spectra: contours 1.2 (dotted), 1.5 (solid); reference shading 0.5-2.75, darker = stronger. (c,d) Kelvin and MJO OLR propagation: reference shading with colour bars in W/m2; thick contours at +/-2 (Kelvin) and +/-3 (MJO), with signed labels. Propagation line styles identify datasets: GC5 solid blue, GC6 dashed red, reference dotted black. Latitude means: 5S-5N Kelvin, 10S-10N MJO; event-centred daily lags, 40-180E. Composite panels remain provisional owing to the source edge-event issue. (e,f) Unchanged area-mean local rainfall/OLR wave-to-total-variance percentages for the SEA wave box.")
    paragraph("Spectral organisation is not uniformly improved.",
              "Both variables retain weaker coherent spectral peaks than the reference. Over the displayed spectral domain, GC6 symmetric power/background RMSE improves 8% for OLR but worsens 3% for rainfall. Seasonal 30-60-day, wavenumber 1-3 east/west power ratios give a clearer MJO penalty: summer rainfall falls 1.86 to 1.74 (reference 2.38), and OLR 1.74 to 1.50 (2.37). Winter rainfall changes little, 1.91 to 1.96 (4.93), while OLR worsens 2.04 to 1.75 (4.55). Background-normalised spectra assess organisation, not absolute amplitude.")
    paragraph("Propagation gains are mainly OLR-specific.",
              "Saved wave composites show eastward Kelvin/MJO and westward Rossby structures. GC6 rainfall-composite RMSE worsens by 3%, 1% and 3% for Kelvin, MJO and Rossby respectively; MJO pattern correlation declines 0.95 to 0.89. OLR composite RMSE improves for MJO (7%) and Rossby (11%), but not Kelvin (+5%). These event-conditioned changes are provisional. Independently, MJO-filtered lag-correlation RMSE improves 9% for OLR but worsens 4% for rainfall. Both simulations lose much of the observed downstream coherence beyond the Indian Ocean.")
    paragraph("OLR fractions do not remove the rainfall limitation.",
              "NDJF MJO OLR contribution increases 6.37% to 7.31%, still well below 12.62% observed; NDJF Kelvin improves 4.96% to 5.97% (7.47%). JJAS Rossby excess decreases 14.05% to 12.02% (7.92%). These gains support improved selected convective fractions, not a general improvement in coupled circulation, rainfall propagation or surge impacts.")
    path = base.HERE / "GC6_focused_MJO_SEA_assessment.docx"
    report.save(path)
    reopened = Document(path)
    assert len(reopened.inline_shapes) == 2
    assert len(reopened.element.xpath("//w:br[@w:type='page']")) == 2
    print(f"PASS: Word report saved with two figures and two detailed-summary pages: {path}")


if __name__ == "__main__":
    import sys

    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9})
    self_test()
    if "--self-test" not in sys.argv:
        if "--document-only" in sys.argv:
            results = json.loads((base.HERE / "focused_evidence.json").read_text())
        elif "--figures-only" in sys.argv:
            results = json.loads((base.HERE / "focused_evidence.json").read_text())
            surge_figure(results)
            wave_figure(results)
        else:
            results = evidence()
            surge_figure(results)
            wave_figure(results)
        document(results)