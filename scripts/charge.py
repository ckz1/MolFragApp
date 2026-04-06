import streamlit as st
import pandas as pd 
import matplotlib.pyplot as plt 
from glob import glob 

import warnings

import numpy as np
from matplotlib.collections import LineCollection


def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)

# csv_file_list = sorted(glob("*let_*/TRAJ_*/*.csv"))
# csv_file = st.radio(
    # "选择数据文件",
    # csv_file_list,
    # horizontal = True   
# )

# df = pd.read_csv(csv_file)

# # st.write(df)

# # df.columns

# x_index = 'time'
# x = df[x_index]
# y_filter = st.text_input("数据关键词", "chg")
# if y_filter == '':
    # y_filter = '-'
# y_indices = st.multiselect(
    # "选择数据",
    # df.columns,
    # [*[item for item in df.columns if y_filter in item]],
# )

# fig, ax = plt.subplots()
# ax.plot(x,df[y_indices])

# ax.set_xlabel(x_index)
# ax.set_ylabel('Charge')
# ax.legend(y_indices)

# st.pyplot(fig)

# csv_file_list = sorted(glob("*let_*/TRAJ_*/*.csv"))
csv_file_list = sorted(glob("*let_*/TRAJ_*/*.csv"))
fig, ax = plt.subplots()
charge_total = 2
for csv_file in csv_file_list:
    df = pd.read_csv(csv_file)
    for i in range(1,2):
        d1, d2, d3 = df[f's{i}-chg1']/charge_total, df[f's{i}-chg2']/charge_total, df[f's{i}-chg3']/charge_total
        x = (d3 - d2) / 3**0.5
        y = d1 
    
        if d3.iloc[-1] < 0.1:
            # ax.plot(x,y,label=f"s{i}",c='gray',alpha=0.6)
            ax.plot(x.iloc[0],y.iloc[0],'r*')
            ax.plot(x.iloc[-1],y.iloc[-1],'bD')
            
            lines = colored_line(x, y, np.linspace(0, 1, x.size), ax, linewidth=1, cmap="plasma")
            
        # print(x.iloc[-1],y.iloc[-1])

ax.plot([1/3**0.5,0,-1/3**0.5,1/3**0.5],[0,1,0,0],'k-')

ax.text(0,0,'Charge #1', ha='center',va='top')
ax.text(1/3**0.5/2,1/2,'Charge #2')
ax.text(-1/3**0.5/2,1/2,'Charge #3',ha='right')

ax.axis('equal')
# ax.legend()
ax.set_xticks([])
ax.set_yticks([])

st.pyplot(fig)