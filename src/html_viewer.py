import os

import jinja2


def _run_plot_dir(out_plot_dir, config):
    """Return the directory used for a run's plots, matching the YAML output layout."""
    if config is None:
        return None

    runid = config.get('runid')
    dir_label = config.get('dir_label', config.get('label', '')).replace(' ', '_')

    options = []
    if runid:
        options.append(os.path.join(out_plot_dir, f'{runid}_{dir_label}'))
    options.extend([
        os.path.join(out_plot_dir, dir_label),
        os.path.join(out_plot_dir, runid) if runid else None,
        out_plot_dir,
    ])

    for path in options:
        if path and os.path.isdir(path):
            return path

    if runid:
        return os.path.join(out_plot_dir, f'{runid}_{dir_label}')
    return os.path.join(out_plot_dir, dir_label)


def generate_html_template(title=None, images=None, html_filename='test.html'):
    """Write a simple HTML gallery showing rows of plots."""
    template = jinja2.Template("""
        <!DOCTYPE html>
        <html>
        <head>
            <title>{{ title }}</title>
            <style>
                :root {
                    --bg: #f4f7fb;
                    --panel: #ffffff;
                    --panel-alt: #eef3fa;
                    --border: #d9e2f0;
                    --text: #1f2a37;
                    --muted: #5e6c82;
                    --shadow: 0 10px 24px rgba(30, 41, 59, 0.08);
                    --accent: #2d6cdf;
                }

                * { box-sizing: border-box; }

                body {
                    margin: 0;
                    padding: 32px 24px 48px;
                    font-family: Arial, Helvetica, sans-serif;
                    background: linear-gradient(180deg, #f8fafc 0%, #eef4fb 100%);
                    color: var(--text);
                }

                h1 {
                    margin: 0 0 24px;
                    font-size: 2rem;
                    font-weight: 700;
                    letter-spacing: -0.03em;
                }

                .container {
                    width: 100%;
                    max-width: 1600px;
                }

                .row {
                    display: flex;
                    flex-wrap: wrap;
                    gap: 20px;
                    margin-bottom: 22px;
                    padding: 18px;
                    background: rgba(255,255,255,0.55);
                    border: 1px solid var(--border);
                    border-radius: 14px;
                    box-shadow: var(--shadow);
                }

                .column {
                    flex: 1 1 300px;
                    min-width: 280px;
                    display: flex;
                    align-items: stretch;
                }

                .image {
                    width: 100%;
                    height: auto;
                    display: block;
                    background: var(--panel);
                    border: 1px solid var(--border);
                    border-radius: 10px;
                    box-shadow: 0 6px 18px rgba(15, 23, 42, 0.08);
                }

                .empty-cell {
                    width: 100%;
                    min-height: 220px;
                    border: 1px dashed #c5d0df;
                    border-radius: 10px;
                    background: repeating-linear-gradient(
                        135deg,
                        var(--panel-alt),
                        var(--panel-alt) 12px,
                        #f7f9fc 12px,
                        #f7f9fc 24px
                    );
                    box-shadow: inset 0 0 0 1px rgba(255,255,255,0.6);
                }
            </style>
        </head>
        <body>
            <h1>{{ title }}</h1>
            <div class="container">
                {% for image_row in images %}
                    <div class="row">
                        {% for image in image_row %}
                            {% if image %}
                                <div class="column">
                                    <img src="{{ image }}" alt="Plot" class="image">
                                </div>
                            {% else %}
                                <div class="column">
                                    <div class="empty-cell" aria-label="No image available"></div>
                                </div>
                            {% endif %}
                        {% endfor %}
                    </div>
                {% endfor %}
            </div>
        </body>
        </html>
    """)

    html_dir = os.path.dirname(os.path.abspath(html_filename))
    relative_images = [
        [os.path.relpath(image, html_dir) if image else None for image in image_row]
        for image_row in images
    ]
    html_content = template.render(title=title, images=relative_images)
    with open(html_filename, 'w') as handle:
        handle.write(html_content)

    print(f'HTML file {html_filename} has been generated.')


def _run_output_label(config):
    """Return the `<runid>_<dir_label>` label used for a run's output folder."""
    runid = str(config.get('runid', '') or '')
    dir_label = str(config.get('dir_label', config.get('label', '')) or '').replace(' ', '_')
    if runid and dir_label:
        return f'{runid}_{dir_label}'
    return runid or dir_label


def _run_prefixes(config):
    """Return the filename prefixes a run may use, longest first."""
    runid = str(config.get('runid', '') or '')
    dir_label = str(config.get('dir_label', config.get('label', '')) or '').replace(' ', '_')

    prefixes = []
    if runid and dir_label:
        prefixes.append(f'{runid}_{dir_label}')
    if dir_label:
        prefixes.append(dir_label)
    if runid:
        prefixes.append(runid)
    return sorted(set(prefixes), key=len, reverse=True)


def _plot_key(stem, prefixes):
    """Strip a run prefix from a plot stem so plots can be matched across runs."""
    for prefix in prefixes:
        if stem.startswith(f'{prefix}_'):
            return stem[len(prefix) + 1:]
    return None


def _plot_index(out_plot_dir, config):
    """Map prefix-stripped plot names to the PNG files available for a run."""
    if config is None:
        return {}

    figure_dir = _run_plot_dir(out_plot_dir, config)
    if not figure_dir or not os.path.isdir(figure_dir):
        return {}

    prefixes = _run_prefixes(config)
    index = {}
    for filename in sorted(os.listdir(figure_dir)):
        if not filename.lower().endswith('.png'):
            continue

        key = _plot_key(os.path.splitext(filename)[0], prefixes)
        if not key or key == 'WH04':
            continue

        path = os.path.abspath(os.path.join(figure_dir, filename))
        existing = index.get(key)
        # Prefer the newest file when old and new naming conventions coexist.
        if existing is None or os.path.getmtime(path) > os.path.getmtime(existing):
            index[key] = path

    return index


def get_obs_basenames(out_plot_dir, config):
    """Return shared plot names for a run, without the runid or label prefix."""
    return sorted(_plot_index(out_plot_dir, config))


def generate_html(out_plot_dir, control=None, expt=None, obs=None):
    """Generate a simple HTML summary page for the available run plots."""
    configs = [config for config in [control, expt, obs] if config is not None]
    if not configs:
        raise ValueError('At least one run config is required to generate the HTML gallery.')

    os.makedirs(out_plot_dir, exist_ok=True)

    reference_config = obs if obs is not None else configs[0]
    indices = [_plot_index(out_plot_dir, config) for config in configs]
    reference_index = indices[configs.index(reference_config)]

    fignames = sorted(reference_index)
    if not fignames:
        print(f'No PNG files found under {out_plot_dir}.')
        return

    images = []
    for figname in fignames:
        row = [index.get(figname) for index in indices]
        images.append(row)

    if not images:
        print('No matching plot rows were found for the configured runs.')
        return

    page_title = ' vs '.join(config['dir_label'] for config in configs)
    if expt is not None and control is not None:
        html_name = f'{_run_output_label(expt)}_v_{_run_output_label(control)}.html'
    elif len(configs) >= 2:
        html_name = f'{_run_output_label(configs[1])}_v_{_run_output_label(configs[0])}.html'
    else:
        html_name = f'{_run_output_label(configs[0])}.html'

    html_file = os.path.join(out_plot_dir, html_name)
    generate_html_template(title=page_title, images=images, html_filename=html_file)
