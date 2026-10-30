import argparse
import os
import logging

import matplotlib
matplotlib.use('Agg', force=True)

import do_mjo
import do_sea
from src import html_viewer
from src import retrieve_SEA_data


SEA_PLOTS_DIR = '/data/users/prince.xavier/SEAPy_climate/SEAPy_output'
MJO_PLOTS_DIR = '/data/users/prince.xavier/SEAPy_climate/MJOPy_output'
RUN_CONFIG_FILE = os.path.join(os.path.dirname(__file__), 'sea_runs.yaml')

os.environ.setdefault('SEAPY_OUTPUT_DIR', SEA_PLOTS_DIR)


logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
logger = logging.getLogger(__name__)


def _parse_value(raw_value):
    value = raw_value.strip()
    if value.startswith("'") and value.endswith("'"):
        return value[1:-1]
    if value.startswith('"') and value.endswith('"'):
        return value[1:-1]
    return value


def load_runs_from_yaml(config_path):
    runs = {}
    current_run_name = None

    logger.info('Loading run configuration from %s', config_path)

    with open(config_path, 'r') as handle:
        for raw_line in handle:
            line = raw_line.rstrip('\n')
            stripped = line.strip()

            if not stripped or stripped.startswith('#'):
                continue

            if stripped == 'runs:':
                continue

            if line.startswith('  ') and not line.startswith('    '):
                if not stripped.endswith(':'):
                    raise ValueError(f'Invalid run section in {config_path}: {line}')
                current_run_name = stripped[:-1]
                runs[current_run_name] = {}
                continue

            if line.startswith('    ') and current_run_name is not None:
                if ':' not in stripped:
                    raise ValueError(f'Invalid run field in {config_path}: {line}')
                key, value = stripped.split(':', 1)
                runs[current_run_name][key.strip()] = _parse_value(value)
                continue

            raise ValueError(f'Unsupported line in {config_path}: {line}')

    return runs


def normalize_run(run):
    normalized = dict(run)
    normalized['label_plot'] = normalized.get('label_plot', normalized['label'])
    normalized['dir_label'] = normalized.get('dir_label', normalized['label'].replace(' ', '_'))
    return normalized


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Run the SEA and MJO diagnostics for any combination of runs.',
        epilog=(
            'Examples:\n'
            '  main.py                                   run everything for all three runs\n'
            '  main.py --runs obs                        only the observations\n'
            '  main.py --runs control experiment         only the two models\n'
            '  main.py --runs obs --sea-levels none      MJO only (no SEA levels selected)\n'
            '  main.py --mjo-levels none                 skip MJO diagnostics\n'
            '  main.py --mjo-levels 1 2 --force          redo MJO levels 1 and 2 for all runs\n'
            '  main.py --html-only                       just rebuild the HTML galleries'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--runs', nargs='+', choices=['control', 'experiment', 'obs', 'all'],
                        default=['all'], help='which runs to process (default: all)')
    parser.add_argument('--mjo-levels', nargs='*', choices=['1', '2', '3', '4', 'none'],
                        default=['1', '2', '3', '4'],
                        help='MJO levels to run; use "none" or pass no values to skip MJO')
    parser.add_argument('--sea-levels', nargs='*', choices=['cs', 'eqw', 'extreme', 'none'],
                        default=['cs', 'eqw', 'extreme'],
                        help='SEA levels to run; use "none" or pass no values to skip SEA')
    parser.add_argument('--force', action='store_true',
                        help='process runs even when all expected plots already exist')
    parser.add_argument('--no-retrieve', action='store_true',
                        help='skip the model data retrieval step')
    parser.add_argument('--html-only', action='store_true',
                        help='only regenerate the HTML galleries from existing plots')
    args = parser.parse_args(argv)
    for option in ('mjo_levels', 'sea_levels'):
        levels = getattr(args, option)
        if 'none' in levels:
            if len(levels) > 1:
                parser.error(f'--{option.replace("_", "-")} none cannot be combined with other levels')
            setattr(args, option, [])
    return args


def select_runs(configs, requested):
    """Return the (name, config) pairs requested on the command line, in a fixed order."""
    if 'all' in requested:
        return list(configs.items())
    return [(name, configs[name]) for name in configs if name in requested]


def _output_label(run):
    return f"{run['runid']}_{run['dir_label'].replace(' ', '_')}"


def runs_needing_work(selected, sea_levels, mjo_levels, force=False):
    """Drop runs whose expected plots already exist, before any data is opened."""
    pending = []
    for name, run in selected:
        if force:
            pending.append((name, run))
            continue

        output_label = _output_label(run)
        missing = do_sea.missing_plot_files(
            run,
            os.path.join(SEA_PLOTS_DIR, output_label),
            os.path.join(MJO_PLOTS_DIR, output_label),
            cs_level1='cs' in sea_levels,
            eqw_level2='eqw' in sea_levels,
            extreme_level3='extreme' in sea_levels,
        )
        missing += do_mjo.missing_plot_files(
            run,
            os.path.join(MJO_PLOTS_DIR, output_label),
            level1='1' in mjo_levels,
            level2='2' in mjo_levels,
            level3='3' in mjo_levels,
            level4_prop='4' in mjo_levels,
        )

        if missing:
            logger.info('%s has %s missing plots; it will be processed', run['runid'], len(missing))
            pending.append((name, run))
        else:
            logger.info('%s already has all expected plots; skipping', run['runid'])

    return pending


def main(argv=None):
    args = parse_args(argv)

    configs = {name: normalize_run(run) for name, run in load_runs_from_yaml(RUN_CONFIG_FILE).items()}
    control = configs.get('control')
    experiment = configs.get('experiment')
    obs = configs.get('obs')

    varnames = ['U850', 'U200', 'OLR', 'V850', 'PRECIP', 'SST']
    sea_levels = set(args.sea_levels)
    mjo_levels = set(args.mjo_levels)

    logger.info('SEA output directory: %s', SEA_PLOTS_DIR)
    logger.info('MJO output directory: %s', MJO_PLOTS_DIR)

    if not args.html_only:
        selected = select_runs(configs, set(args.runs))
        if not selected:
            logger.warning('No runs selected; nothing to do.')
            return

        logger.info('Runs selected: %s', ', '.join(name for name, _ in selected))
        logger.info('MJO levels: %s', ', '.join(sorted(mjo_levels)) or 'none')
        logger.info('SEA levels: %s', ', '.join(sorted(sea_levels)) or 'none')

        logger.info('Checking for existing plots before retrieving data')
        pending = dict(runs_needing_work(selected, sea_levels, mjo_levels, force=args.force))

        if not pending:
            logger.info('All expected plots are present; skipping data retrieval and computation')
        else:
            pending_control = pending.get('control')
            pending_experiment = pending.get('experiment')
            pending_obs = pending.get('obs')

            if not args.no_retrieve and (pending_control or pending_experiment):
                logger.info('Ensuring required SEA data are present before computation')
                retrieve_SEA_data.model_data_retrieve(
                    varnames,
                    control=pending_control,
                    expt=pending_experiment,
                    netcdf=True,
                )

            if mjo_levels:
                # MJO runs first because the SEA extremes diagnostics use its RMM index files.
                do_mjo.mjo_compute(
                    control=pending_control,
                    expt=pending_experiment,
                    obs=pending_obs,
                    level1='1' in mjo_levels,
                    level2='2' in mjo_levels,
                    level3='3' in mjo_levels,
                    level4_prop='4' in mjo_levels,
                    mjo_plots_dir=MJO_PLOTS_DIR,
                )
                logger.info('Completed MJO workflow')

            if sea_levels:
                do_sea.sea_compute(
                    varnames,
                    control=pending_control,
                    expt=pending_experiment,
                    obs=pending_obs,
                    cs_level1='cs' in sea_levels,
                    eqw_level2='eqw' in sea_levels,
                    extreme_level3='extreme' in sea_levels,
                    sea_plots_dir=SEA_PLOTS_DIR,
                    mjo_plots_dir=MJO_PLOTS_DIR,
                )
                logger.info('Completed SEA workflow')

    # The galleries always use every configured run so comparison pages stay complete.
    logger.info('Generating MJO HTML viewer in %s', MJO_PLOTS_DIR)
    html_viewer.generate_html(MJO_PLOTS_DIR, control=control, expt=experiment, obs=obs)
    logger.info('Generating SEA HTML viewer in %s', SEA_PLOTS_DIR)
    html_viewer.generate_html(SEA_PLOTS_DIR, control=control, expt=experiment, obs=obs)


if __name__ == '__main__':
    main()
