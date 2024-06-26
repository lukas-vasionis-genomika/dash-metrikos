import numpy as np
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
import seaborn as sns
import pandas as pd
from .utils import nullity_filter, nullity_sort
import warnings


def matrix(
    df, filter=None, n=0, p=0, sort=None, figsize=(25, 10), width_ratios=(15, 1),
    color=(0.25, 0.25, 0.25), fontsize=16, labels=None, label_rotation=45, sparkline=True,
    freq=None, ax=None, row_labels=None,ax_title=None
):
    """
    A matrix visualization of the nullity of the given DataFrame.

    :param df: The `DataFrame` being mapped.
    :param filter: The filter to apply to the heatmap. Should be one of "top", "bottom", or None (default).
    :param n: The max number of columns to include in the filtered DataFrame.
    :param p: The max percentage fill of the columns in the filtered DataFrame.
    :param sort: The row sort order to apply. Can be "ascending", "descending", or None.
    :param figsize: The size of the figure to display.
    :param fontsize: The figure's font size. Default to 16.
    :param labels: Whether or not to display the column names. Defaults to the underlying data labels when there are
        50 columns or less, and no labels when there are more than 50 columns.
    :param label_rotation: What angle to rotate the text labels to. Defaults to 45 degrees.
    :param sparkline: Whether or not to display the sparkline. Defaults to True.
    :param width_ratios: The ratio of the width of the matrix to the width of the sparkline. Defaults to `(15, 1)`.
        Does nothing if `sparkline=False`.
    :param color: The color of the filled columns. Default is `(0.25, 0.25, 0.25)`.
    :return: The plot axis.
    """
    df = nullity_filter(df, filter=filter, n=n, p=p)
    df = nullity_sort(df, sort=sort, axis='columns')

    height = df.shape[0]
    width = df.shape[1]

    # z is the color-mask array, g is a NxNx3 matrix. Apply the z color-mask to set the RGB of each pixel.
    z = df.notnull().values
    g = np.zeros((height, width, 3), dtype=np.float32)

    g[z < 0.5] = [1, 1, 1]
    g[z > 0.5] = color

    # Set up the matplotlib grid layout. A unary subplot if no sparkline, a left-right splot if yes sparkline.
    if ax is None:
        plt.figure(figsize=figsize)
        if sparkline:
            gs = gridspec.GridSpec(1, 2, width_ratios=width_ratios)
            gs.update(wspace=0.08)
            ax1 = plt.subplot(gs[1])
        else:
            gs = gridspec.GridSpec(1, 1)
        ax0 = plt.subplot(gs[0])
    else:
        if sparkline is not False:
            warnings.warn(
                "Plotting a sparkline on an existing axis is not currently supported. "
                "To remove this warning, set sparkline=False."
            )
            sparkline = False
        ax0 = ax

    # Create the nullity plot.
    ax0.imshow(g, interpolation='none')

    # Remove extraneous default visual elements.
    ax0.set_aspect('auto')
    ax0.grid(visible=False)
    ax0.xaxis.tick_top()
    ax0.xaxis.set_ticks_position('none')
    ax0.yaxis.set_ticks_position('none')
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['bottom'].set_visible(False)
    ax0.spines['left'].set_visible(False)

    # Set up and rotate the column ticks. The labels argument is set to None by default. If the user specifies it in
    # the argument, respect that specification. Otherwise display for <= 50 columns and do not display for > 50.
    if type(labels)==list:
        ha = 'left'
        ax0.set_xticks(list(range(0, width)))
        ax0.set_xticklabels(labels, rotation=label_rotation, ha=ha, fontsize=fontsize)

    # elif labels or (labels is None and len(df.columns) <= 50):
    #     ha = 'left'
    #     ax0.set_xticks(list(range(0, width)))
    #     ax0.set_xticklabels(list(df.columns), rotation=label_rotation, ha=ha, fontsize=fontsize)
    else:
        ax0.set_xticks([])

    # Adds Timestamps ticks if freq is not None, else set up the two top-bottom row ticks.
    # if freq:
    #     ts_list = []
    #
    #     if type(df.index) == pd.PeriodIndex:
    #         ts_array = pd.date_range(df.index.to_timestamp().date[0],
    #                                  df.index.to_timestamp().date[-1],
    #                                  freq=freq).values
    #
    #         ts_ticks = pd.date_range(df.index.to_timestamp().date[0],
    #                                  df.index.to_timestamp().date[-1],
    #                                  freq=freq).map(lambda t:
    #                                                 t.strftime('%Y-%m-%d'))
    #
    #     elif type(df.index) == pd.DatetimeIndex:
    #         ts_array = pd.date_range(df.index[0], df.index[-1],
    #                                  freq=freq).values
    #
    #         ts_ticks = pd.date_range(df.index[0], df.index[-1],
    #                                  freq=freq).map(lambda t:
    #                                                 t.strftime('%Y-%m-%d'))
    #     else:
    #         raise KeyError('Dataframe index must be PeriodIndex or DatetimeIndex.')
    #     try:
    #         for value in ts_array:
    #             ts_list.append(df.index.get_loc(value))
    #     except KeyError:
    #         raise KeyError('Could not divide time index into desired frequency.')
    #
    #     ax0.set_yticks(ts_list)
    #     ax0.set_yticklabels(ts_ticks, fontsize=int(fontsize / 16 * 20), rotation=0)
    if row_labels is list:
        ax0.set_yticks(row_labels)
        ax0.set_yticklabels(row_labels, fontsize=fontsize, rotation=0)

    # Create the inter-column vertical grid.
    in_between_point = [x + 0.5 for x in range(0, width - 1)]
    for in_between_point in in_between_point:
        ax0.axvline(in_between_point, linestyle='-', color='white')

    ax0.set_title(ax_title)
    # if sparkline:
    #     # Calculate row-wise completeness for the sparkline.
    #     completeness_srs = df.notnull().astype(bool).sum(axis=1)
    #     x_domain = list(range(0, height))
    #     y_range = list(reversed(completeness_srs.values))
    #     min_completeness = min(y_range)
    #     max_completeness = max(y_range)
    #     min_completeness_index = y_range.index(min_completeness)
    #     max_completeness_index = y_range.index(max_completeness)
    #
    #     # Set up the sparkline, remove the border element.
    #     ax1.grid(visible=False)
    #     ax1.set_aspect('auto')
    #     # GH 25
    #     if int(mpl.__version__[0]) <= 1:
    #         ax1.set_axis_bgcolor((1, 1, 1))
    #     else:
    #         ax1.set_facecolor((1, 1, 1))
    #     ax1.spines['top'].set_visible(False)
    #     ax1.spines['right'].set_visible(False)
    #     ax1.spines['bottom'].set_visible(False)
    #     ax1.spines['left'].set_visible(False)
    #     ax1.set_ymargin(0)
    #
    #     # Plot sparkline---plot is sideways so the x and y axis are reversed.
    #     ax1.plot(y_range, x_domain, color=color)
    #
    #     if labels:
    #         # Figure out what case to display the label in: mixed, upper, lower.
    #         label = 'Data Completeness'
    #         if str(df.columns[0]).islower():
    #             label = label.lower()
    #         if str(df.columns[0]).isupper():
    #             label = label.upper()
    #
    #         # Set up and rotate the sparkline label.
    #         ha = 'left'
    #         ax1.set_xticks([min_completeness + (max_completeness - min_completeness) / 2])
    #         ax1.set_xticklabels([label], rotation=label_rotation, ha=ha, fontsize=fontsize)
    #         ax1.xaxis.tick_top()
    #         ax1.set_yticks([])
    #     else:
    #         ax1.set_xticks([])
    #         ax1.set_yticks([])
    #
    #     # Add maximum and minimum labels, circles.
    #     ax1.annotate(max_completeness,
    #                  xy=(max_completeness, max_completeness_index),
    #                  xytext=(max_completeness + 2, max_completeness_index),
    #                  fontsize=int(fontsize / 16 * 14),
    #                  va='center',
    #                  ha='left')
    #     ax1.annotate(min_completeness,
    #                  xy=(min_completeness, min_completeness_index),
    #                  xytext=(min_completeness - 2, min_completeness_index),
    #                  fontsize=int(fontsize / 16 * 14),
    #                  va='center',
    #                  ha='right')
    #
    #     ax1.set_xlim([min_completeness - 2, max_completeness + 2])  # Otherwise the circles are cut off.
    #     ax1.plot([min_completeness], [min_completeness_index], '.', color=color, markersize=10.0)
    #     ax1.plot([max_completeness], [max_completeness_index], '.', color=color, markersize=10.0)
    #
    #     # Remove tick mark (only works after plotting).
    #     ax1.xaxis.set_ticks_position('none')
    return ax0