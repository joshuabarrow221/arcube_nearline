#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

DATES = ["04/09", "05/09", "06/09", "07/09", "08/09", "09/09", "10/09"] # to make it more precise use time of measurement
PURITY_MONITOR_LIFETIMES_us = [114, 66.5, 156, 182, 196, 206, 215]
# add errors
GAS_ANALYSER_O2 = [30.7, 17.6, 15.0, 15.0, 14.9, 14.9, 159]
GAS_ANALYSER_H2O = [69.3, 60.7, 54.7, 50.8, 48.5, 47.0, 219]

# Check how correct this is correct ^^'
def lifetime_gas_analyser_us(O2, H2O):
    lifetime_O2 = 0.299 / O2
    lifetime_H2O = 17 / H2O
    one_over_tau = (1 / lifetime_O2) + (1 / lifetime_H2O)
    return (1 / one_over_tau) * 1000

def plot_lifetime(lifetime_monitor, lifetime_gas, dates, projection_days=15):
    x = np.array([datetime.strptime(f"2026/{date}", "%Y/%d/%m") for date in dates])
    y_monitor = np.array(lifetime_monitor)
    y_gas = np.array(lifetime_gas)
    anode_start_index = 2
    x_days = np.array([(d - x[0]).total_seconds() / 86400 for d in x])
    fit_monitor = np.polyfit(x_days[anode_start_index:], y_monitor[anode_start_index:], 1)
    fit_gas = np.polyfit(x_days[anode_start_index:], y_gas[anode_start_index:], 1)
    future_x = np.array([x[-1] + timedelta(days=i) for i in range(1, projection_days + 1)])
    future_days = np.array([(d - x[0]).total_seconds() / 86400 for d in future_x])
    future_monitor = np.polyval(fit_monitor, future_days)
    future_gas = np.polyval(fit_gas, future_days)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, y_monitor, "o-", color="tab:blue", linewidth=2, markersize=7, label="Purity monitor")
    ax.plot(x, y_gas, "o-", color="tab:orange", linewidth=2, markersize=7, label="Gas analyser")
    ax.plot(future_x, future_monitor, "--", color="tab:blue", linewidth=2, alpha=0.7, label="Monitor projection")
    ax.plot(future_x, future_gas, "--", color="tab:orange", linewidth=2, alpha=0.7, label="Gas analyser projection")
    dashed_date = x[0] + (x[1] - x[0]) / 2
    ax.axvline(dashed_date, color="black", linestyle="--", linewidth=1.5)
    ymax = max(np.max(y_monitor), np.max(y_gas), np.max(future_monitor), np.max(future_gas))
    ax.text(dashed_date, ymax * 0.95, "Anode signal starts", rotation=90, verticalalignment="top", horizontalalignment="right", fontsize=11)
    ax.set_xlabel("Date", fontsize=12)
    ax.set_ylabel("Lifetime [µs]", fontsize=12)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.autofmt_xdate()
    plt.tight_layout()
    plt.savefig("lifetime_history.png")


def main():
    lifetime_gas = lifetime_gas_analyser_us(np.array(GAS_ANALYSER_O2), np.array(GAS_ANALYSER_H2O))
    print("Gas analyser lifetimes [us]:")
    print(lifetime_gas)
    plot_lifetime(PURITY_MONITOR_LIFETIMES_us, lifetime_gas, DATES, projection_days=20)


if __name__ == '__main__':
    main()
