import os

import matplotlib.pyplot as plt


def save_figure(filename, *args, **kwargs):
    plt.savefig(filename, *args, **kwargs)
    root, extension = os.path.splitext(filename)
    if extension.lower() == '.png':
        plt.savefig(f'{root}.pdf', *args, **kwargs)