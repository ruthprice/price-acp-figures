# time_series_pdf.py
# Plots the times series and PDFs of aerosol during AO2018
# eg for figure 2
# =====================================================================
import matplotlib.pyplot as plt
import figure_props as fprops
import numpy as np
import matplotlib.lines as mlines
import datetime as dt
import matplotlib.dates as mdates
# =====================================================================
def plot_time_series_pdfs(suites, model_N, model_T, obs_N, obs_T, obs_stdev,
                          hist_obs, hist_model, pdf_bins_mid, filename):
    sizes = ['2.5--15 nm', '15--100 nm', '100--500 nm']
    n_Dp_bins = len(sizes)
    fig = plt.figure(figsize=(16*fprops.cm,14*fprops.cm), dpi=fprops.dpi)
    widths = [2, 1]
    heights = [1, 1, 1]
    spec = fig.add_gridspec(ncols=2, nrows=3,
                            width_ratios=widths,
                            height_ratios=heights,
                            hspace=0.4, wspace=0.125)
    lower = [obs_N[n] - 0.5*obs_stdev[n] for n in np.arange(n_Dp_bins)]
    upper = [obs_N[n] + 0.5*obs_stdev[n] for n in np.arange(n_Dp_bins)]
    for i in np.arange(n_Dp_bins):
        # PLOT TIMES SERIES
        ax1 = fig.add_subplot(spec[i,0])
        legend_lines = []
        # plot obs
        ax1.fill_between(obs_T[i], lower[i], upper[i], color='lightgrey', alpha=0.75)
        ax1.plot(obs_T[i], obs_N[i], color='black', linewidth=fprops.thick_line)
        legend_lines.append(mlines.Line2D([],[],linestyle='solid',color='black',label='Obs'))
        # plot model output
        for s,suite in enumerate(suites):
            colour = fprops.colours[suite]
            linewidth = fprops.linewidths[suite]
            label = fprops.suite_labels[suite]
            ax1.plot(model_T[suite], model_N[suite][:,i],
                     color=colour, linewidth=linewidth)
            legend_lines.append(mlines.Line2D([],[], color=colour, label=label, linewidth=linewidth))
        # mark NPF events
        for j,event in enumerate(fprops.ao2018_npf_events):
            ax1.hlines(fprops.npf_event_y[i],event[0],event[1],
                       linewidth=0.8, color='red', linestyle=(0,(1.8,0.9)))
            ax1.plot(event[0], fprops.npf_event_y[i], marker='|', ms=2.5, color='red')
            ax1.plot(event[1], fprops.npf_event_y[i], marker='|', ms=2.5,color='red')
        # make pretty axes
        ax1.tick_params(axis='y', labelsize=fprops.ax_fs)
        ax1.tick_params(axis='y', which='major', pad=1)
        ax1.set_yscale('log')
        if i == 0:
            ax1.set_ylim(bottom=1e-6,top=5e4)
        if i == 1:
            ax1.set_ylabel('Particle concentration [cm$^{-3}$]',
                           fontsize=fprops.ax_label_fs, labelpad=3)
        if i == 2:
            ax1.set_xlabel('Day of year', fontsize=fprops.ax_label_fs)
        left_lim_date = dt.datetime(2018, 8, 2)
        right_lim_date = dt.datetime(2018, 9, 20)
        ax1.set_xlim(left=left_lim_date, right=right_lim_date)
        ax1.tick_params(axis='x', labelsize=fprops.ax_fs)
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%d-%m'))
        # make pretty legend
        if i == 0:
            ax1.legend(handles=legend_lines, fontsize=fprops.legend_fs,
                       loc='lower center', bbox_to_anchor=(0.5,0),
                       ncol=3, columnspacing=1, borderpad=0.2,
                       handletextpad=0.3, labelspacing=0.5,
                       handlelength=1.5)
        # make plot titles
        ax1.set_title(sizes[i], fontsize=fprops.label_fs)
        subfig_labels_text = ['(a)','(d)','(g)']
        ax1.set_title(subfig_labels_text[i], loc='left',
                      fontsize=fprops.label_fs, y=1)
        plt.grid()
        # PLOT PDFs
        subfig_labels_text = np.array([['(b)','(c)'],['(e)','(f)'],['(h)','(i)']])
        spec_pdf = spec[i,1].subgridspec(ncols=2, nrows=1, wspace=0.125)
        for p,period in enumerate(['Melt','Freeze']):
            ax = fig.add_subplot(spec_pdf[p])
            # plot obs
            ax.plot(pdf_bins_mid[i], hist_obs[i,p], color='k', label='Obs', linewidth=fprops.thick_line)
            # plot model output
            for s,suite in enumerate(suites):
                colour = fprops.colours[suite]
                linewidth = fprops.linewidths[suite]
                ax.plot(pdf_bins_mid[i], hist_model[suite][i,p],
                        color=colour, linewidth=linewidth)
            ax.grid()
            # make pretty axes
            ax.tick_params(axis='both', which='major', labelsize=fprops.ax_fs)
            ax.tick_params(axis='both', which='minor', labelsize=fprops.ax_fs)
            ax.tick_params(axis='y', which='major', pad=1)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(1e-1,1e4)
            ax.set_ylim(1e-6, 1e0)
            if i == 2:
                ax.set_xlabel("N [cm$^{-3}$]", fontsize=fprops.ax_label_fs)
            if p == 1:
                ax.yaxis.set_ticklabels([])
            # make titles
            ax.text(0.5, 0.95, period, fontsize=fprops.label_fs,
                    va='top', ha='center', transform=ax.transAxes)
            ax.set_title(subfig_labels_text[i,p], loc='left',
                         fontsize=fprops.label_fs, y=1)

    plt.savefig(filename, bbox_inches="tight", facecolor='white', format='pdf')
    plt.close()

def plot_IA_time_series_pdf(obs_IA, obs_T, model_IA, model_T,
                            obs_hist, obs_bins, model_hist, model_bins,
                            filename):
    npf_event_marker_y = 8e6

    fig = plt.figure(figsize=(16*fprops.cm,3*fprops.cm), dpi=300)
    widths = [4,1]
    spec = fig.add_gridspec(ncols=2, nrows=1,
                            width_ratios=widths,
                            hspace=0.4, wspace=0.1)
    ax1 = fig.add_subplot(spec[0])
    # obs
    ax1.plot(obs_T, obs_IA, color='black', label="Observations (3hr mean)", linewidth=fprops.thick_line)
    # model
    ax1.plot(model_T, model_IA, label=fprops.suite_labels['u-cm612'],
             color=fprops.colours['u-cm612'], linewidth=fprops.linewidths['u-cm612'])
    # mark NPF events
    for j,event in enumerate(fprops.ao2018_npf_events):
        ax1.hlines(npf_event_marker_y,event[0],event[1],
                   linewidth=0.8, color='red', linestyle=(0,(1.8,0.9)))
        ax1.plot(event[0], npf_event_marker_y, marker='|', ms=2.5, color='red')
        ax1.plot(event[1], npf_event_marker_y, marker='|', ms=2.5,color='red')
    # make pretty axes
    ax1.set_xlabel('DoY', fontsize=fprops.ax_fs)
    ax1.set_xlim(right=dt.datetime(2018,9,20))
    ax1.set_ylabel('Surface IA conc [cm$^{-3}$]', fontsize=fprops.ax_fs)
    ax1.tick_params(axis='y', which='major', labelsize=fprops.ax_fs)
    ax1.tick_params(axis='y', which='minor', labelsize=fprops.ax_fs)
    ax1.tick_params(axis='y', which='major', pad=1)
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # use DOY on x axis
    left_lim_date = dt.datetime(2018, 8, 1)
    right_lim_date = dt.datetime(2018, 9, 19)
    ax1.set_xlim(left=left_lim_date, right=right_lim_date)
    left,right = plt.xlim()
    doy_left = left_lim_date.timetuple().tm_yday
    doy_right = right_lim_date.timetuple().tm_yday
    x_axis_array = np.arange(left,right+1)
    doy_array = np.arange(doy_left, doy_right+1)
    ax1.set_xticks(x_axis_array[::3])
    ax1.set_xticklabels(doy_array[::3])
    ax1.tick_params(axis='x', which='major', labelsize=fprops.ax_fs)
    ax1.tick_params(axis='x', which='minor', labelsize=fprops.ax_fs)
    ax1.legend(fontsize=fprops.legend_fs)
    ax1.set_title('(a)', loc='left', fontsize=fprops.label_fs, y=1.1)
    # PDFs
    ax2 = fig.add_subplot(spec[1])
    ax2.plot(obs_bins, obs_hist, color='k',
             label='Observations', linewidth=fprops.thick_line)
    ax2.plot(model_bins, model_hist, label=fprops.suite_labels['u-cm612'],
             color=fprops.colours['u-cm612'], linewidth=fprops.linewidths['u-cm612'])
    # Pretty axes
    ax2.set_xlabel("Surface IA conc [cm$^{-3}$]", labelpad=12, fontsize=fprops.ax_fs)
    ax2.tick_params(axis='both', which='major', labelsize=fprops.ax_fs)
    ax2.tick_params(axis='both', which='minor', labelsize=fprops.ax_fs)
    ax2.tick_params(axis='y', which='major', pad=1)
    ax2.set_title('Freeze season', loc='right', fontsize=fprops.label_fs, y=1.1)
    ax2.set_title('(b)', loc='left', fontsize=fprops.label_fs, y=1.1)
    ax1.grid()
    ax2.grid()
    plt.savefig(filename+".pdf", bbox_inches='tight', facecolor='white', format='pdf')
