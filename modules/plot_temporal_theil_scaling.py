# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 2023

@author: Felipe Abril Bermúdez
"""

# Libraries ----
import re
import sys
import warnings
import matplotlib
import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
import matplotlib.colors as mcolors
import matplotlib.transforms as mtransforms

from scipy.optimize import curve_fit

# Global options ----
warnings.filterwarnings("ignore")
pd.options.mode.chained_assignment = None
pd.set_option('display.max_columns', None)

colors = ["blue", "green", "red", "purple", "orange", "brown", "pink", "olive", "gray", "cyan"]

# Estimate tfs parameters (TTS) ----
def temporal_theil_scaling(
    theil_values,
    mean_values,
    tts_coefficient,
    tts_exponent
):
    """Estimation of temporal Theil scaling
    Temporal Theil scaling law:
        theil_values: Theil values from diffusive trajectory time series
        mean_values: Average values from diffusive trajectory time series
    """
    
    mean_values = 1- mean_values / np.max(mean_values)
    tts_value = tts_coefficient * np.power(mean_values, tts_exponent)

    return tts_value

# Plot evolution with TTS parameters ----
def plot_tts_evolution(
    df_tts,
    symbols,
    width,
    height,
    markersize = 2,
    fontsize_labels = 13.5,
    fontsize_legend = 11.5,
    usetex = False,
    n_cols = 4,
    n_x_breaks = 10,
    n_y_breaks = 10,
    fancy_legend = True,
    dpi = 150,
    save_figures = True,
    output_path = "../output_files",
    information_name = "",
    input_generation_date = "2023-03-28"
):
    """Preparation of data for plotting
    Join original data with optimal window size data:
        df_tts: Dataframe with multiple diffusive trajectory time series
        symbols: Symbols of the financial time series plotted
        width: Width of final plot
        height: Height of final plot
        markersize: Marker size as in plt.plot()
        fontsize_labels: Font size in axis labels
        fontsize_legend: Font size in legend
        usetex: Use LaTeX for renderized plots
        n_cols: Number of columns in legend
        n_x_breaks: Number of divisions in x-axis
        n_y_breaks: Number of divisions in y-axis
        fancy_legend: Fancy legend output
        dpi: Dot per inch for output plot
        save_figures: Save figures flag
        output_path: Output path where figures is saved
        information_name: Name of the output plot
        input_generation_date: Date of generation (control version)
    """
    
    # Plot data and define loop over symbols ----
    df_graph = df_tts[df_tts["time_series"].isin(symbols)]
    loop_index = sorted(df_graph["time_series"].unique().tolist())
    
    # Begin plot inputs ----
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "text.usetex": usetex,
            "pgf.rcfonts": False
        }
    )
    
    fig1, ax1 = plt.subplots(len(loop_index), 3)
    fig1.set_size_inches(w = width, h = height)
    counter = 0
        
    for i in loop_index:
        counter_i = 0
        for j in sorted(df_graph[df_graph["time_series"] == i]["sub_time_series"].unique().tolist()):
            # Filter information ----
            df_aux = df_graph[((df_graph["time_series"] == i) & (df_graph["sub_time_series"] == j))]
        
            # Parameters ----
            ave_tts_j = df_aux["mean_value"]
            theil_tts_j = df_aux["norm_sum_theil_index"] / np.max(df_aux["norm_sum_theil_index"])
            rsquared_j = df_aux["rsquared"].unique()[0]
            exponent_tts_j = df_aux["tts_exponent"].unique()[0]
            coefficient_tts_j = df_aux["tts_coefficient"].unique()[0]
            exponent_tts_error_j = df_aux["tts_exponent_error"].unique()[0]
            coefficient_tts_error_j = df_aux["tts_coefficient_error"].unique()[0]

            # Extract parameters of error ----
            tts_mean = temporal_theil_scaling(
                theil_values = theil_tts_j,
                mean_values = ave_tts_j,
                tts_coefficient = coefficient_tts_j,
                tts_exponent = exponent_tts_j
            )
            tts_lower = temporal_theil_scaling(
                theil_values = theil_tts_j,
                mean_values = ave_tts_j,
                tts_coefficient = coefficient_tts_j - coefficient_tts_error_j,
                tts_exponent = exponent_tts_j - exponent_tts_error_j
            )
            tts_upper = temporal_theil_scaling(
                theil_values = theil_tts_j,
                mean_values = ave_tts_j,
                tts_coefficient = coefficient_tts_j + coefficient_tts_error_j,
                tts_exponent = exponent_tts_j + exponent_tts_error_j
            )
            
            tts_mean = tts_mean / np.max(tts_mean)
            tts_lower = tts_lower / np.max(tts_lower)
            tts_upper = tts_upper / np.max(tts_upper)

            # Plot graphs ----
            if len(loop_index) == 1:
                # Plot graph (TTS) ----
                plot_1 = ax1[counter_i].plot(
                    ave_tts_j,
                    theil_tts_j,
                    alpha = 1,
                    zorder = 2,
                    color = "black",
                    marker = "o",
                    linestyle = "",
                    label = "empirical data",
                    markersize = markersize
                )
                ax1[counter_i].plot(ave_tts_j, tts_mean, alpha = 1, zorder = 1, color = colors[counter_i], linewidth = 3, label = "fitting")
                ax1[counter_i].fill_between(
                    ave_tts_j,
                    tts_lower,
                    tts_upper,
                    where = ((tts_upper >= tts_lower) & (tts_upper >= tts_mean) & (tts_mean >= tts_lower)),
                    alpha = 0.19,
                    facecolor = colors[counter_i],
                    interpolate = True
                )
                ax1[counter_i].tick_params(which = "major", direction = "in", top = True, right = True, labelsize = fontsize_labels, length = 12)
                ax1[counter_i].tick_params(which = "minor", direction = "in", top = True, right = True, labelsize = fontsize_labels, length = 6)
                ax1[counter_i].xaxis.set_major_locator(mtick.MaxNLocator(n_x_breaks))
                ax1[counter_i].xaxis.set_minor_locator(mtick.MaxNLocator(4 * n_x_breaks))
                ax1[counter_i].yaxis.set_major_locator(mtick.MaxNLocator(n_y_breaks))
                ax1[counter_i].yaxis.set_minor_locator(mtick.MaxNLocator(5 * n_y_breaks))
                ax1[counter_i].xaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
                ax1[counter_i].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
                ax1[counter_i].tick_params(axis = "x", labelrotation = 90)
                ax1[counter_i].set_xlabel("Mean - {}".format(j.capitalize()), fontsize = fontsize_labels)
                ax1[counter_i].set_ylabel("Normalized Shannon index - {}".format(j.capitalize()), fontsize = fontsize_labels)
                ax1[counter_i].legend(fancybox = fancy_legend, shadow = True, ncol = n_cols, fontsize = fontsize_legend)
                ax1[counter_i].set_title(
                    r"({}) $R^2={}\%$".format(chr(counter_i + 65), round(rsquared_j * 100, 2)),
                    loc = "left",
                    y = 1.005,
                    fontsize = fontsize_labels
                )

            else:
                # Plot graph (TTS) ----
                plot_1 = ax1[counter, counter_i].plot(
                    ave_tts_j,
                    theil_tts_j,
                    alpha = 1,
                    zorder = 2,
                    color = "black",
                    marker = "o",
                    linestyle = "",
                    label = "empirical data",
                    markersize = markersize
                )
                ax1[counter, counter_i].plot(ave_tts_j, tts_mean, alpha = 1, zorder = 1, color = colors[counter_i], linewidth = 3, label = "fitting")
                ax1[counter, counter_i].fill_between(
                    ave_tts_j,
                    tts_lower,
                    tts_upper,
                    where = ((tts_upper >= tts_lower) & (tts_upper >= tts_mean) & (tts_mean >= tts_lower)),
                    alpha = 0.19,
                    facecolor = colors[counter_i],
                    interpolate = True
                )
                ax1[counter, counter_i].tick_params(which = "major", direction = "in", top = True, right = True, labelsize = fontsize_labels, length = 12)
                ax1[counter, counter_i].tick_params(which = "minor", direction = "in", top = True, right = True, labelsize = fontsize_labels, length = 6)
                ax1[counter, counter_i].xaxis.set_major_locator(mtick.MaxNLocator(n_x_breaks))
                ax1[counter, counter_i].xaxis.set_minor_locator(mtick.MaxNLocator(4 * n_x_breaks))
                ax1[counter, counter_i].yaxis.set_major_locator(mtick.MaxNLocator(n_y_breaks))
                ax1[counter, counter_i].yaxis.set_minor_locator(mtick.MaxNLocator(5 * n_y_breaks))
                ax1[counter, counter_i].xaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
                ax1[counter, counter_i].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
                ax1[counter, counter_i].tick_params(axis = "x", labelrotation = 90)
                ax1[counter, counter_i].set_xlabel("Mean - {}".format(j.capitalize()), fontsize = fontsize_labels)
                ax1[counter, counter_i].set_ylabel("Normalized Shannon index - {}".format(j.capitalize()), fontsize = fontsize_labels)
                ax1[counter, counter_i].legend(fancybox = fancy_legend, shadow = True, ncol = n_cols, fontsize = fontsize_legend)
                ax1[counter, counter_i].set_title(
                    r"({}) $R^2={}\%$".format(chr(counter_i + 65), round(rsquared_j * 100, 2)),
                    loc = "left",
                    y = 1.005,
                    fontsize = fontsize_labels
                )

            # Function development ----
            counter_i += 1
            print("Generated plot for {} and time series {}".format(i, j))
        
        counter += 1
    
    fig1.tight_layout()
    if save_figures:
        plt.show()
        fig1.savefig(
            "{}/{}_tts_exponent_evolution_{}.png".format(output_path, information_name, re.sub("-", "", input_generation_date)),
            bbox_inches = "tight",
            facecolor = fig1.get_facecolor(),
            transparent = False,
            pad_inches = 0.03,
            dpi = dpi
        )
    plt.close()
    
    return df_graph
