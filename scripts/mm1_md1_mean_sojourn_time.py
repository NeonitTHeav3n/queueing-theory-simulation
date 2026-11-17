"""Сравнение среднего времени пребывания E[T] в системах M/M/1 и M/D/1.

Проводится:
  1. сравнение "быстрой" симуляции M/M/1 и алгоритма 3.1 (GI/G/K) с теорией
     (Бочаров, Печинкин, «Теория массового обслуживания», стр. 110–114);
  2. сравнение систем M/M/1 и M/D/1;
  3. анализ относительной ошибки по диапазонам загрузки rho;
  4. отношение E[T]_MD1 / E[T]_MM1 (теория и симуляция).

Метрики качества: RMSE, MAE. Графики сохраняются в папку figures/.
"""

import os
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_squared_error, mean_absolute_error

FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams["figure.figsize"] = (9, 5)
sns.set(style="whitegrid")
np.random.seed(42)
random.seed(42)

# Книга Бочаров, Печинкин «Теория массового обслуживания», стр. 110–114:
#   - алгоритм 3.1 моделирования системы GI/G/K;
#   - формулы для E[N] (стр. 114), по которым считается и E[T].
mu = 1.0
num_samples = 200
simulation_time = 250
rho_values = np.linspace(0.01, 0.95, num_samples)


# Алгоритм 3.1 — моделирование системы GI/G/K (Бочаров, Печинкин, стр. 110–114)
def simulate_gi_g_k(F_A, F_B, lam, mean_b, K, N_min=0, N_max=100000, T_max=0, eps=1e-12):
    n = 0
    x_a = 0.0
    x_b = [float("inf")]
    N_a = 0
    T_s = 0.0
    T_n = 0.0

    while N_a < N_max or T_s < T_max:
        T_x = min(x_a, min(x_b))

        x_a = x_a - T_x
        x_b = [xb - T_x for xb in x_b]

        if N_a >= N_min:
            T_n = T_n + n * T_x
            T_s = T_s + T_x

        if x_a <= eps:
            N_a = N_a + 1
            n = n + 1
            x_a = F_A()

            if n <= K:
                x_b.append(F_B())

        elif min(x_b) <= eps:
            j = x_b.index(min(x_b))
            x_b.pop(j)
            n = n - 1

            if n >= K:
                x_b.append(F_B())

        if N_a >= N_max and T_s >= T_max:
            break

    if T_s > 0:
        E_N = T_n / T_s
        E_T = E_N / lam
        E_W = E_T - mean_b
        E_Q = lam * E_W
    else:
        E_N = E_T = E_W = E_Q = 0.0

    return {
        "E[N]": E_N,
        "E[T]": E_T,
        "E[W]": E_W,
        "E[Q]": E_Q
    }


# "Быстрая" симуляция M/M/1 (построение траектории на numpy)
def simulate_mm1_fast(rho, mu=1.0, sim_time=250):
    if rho >= 0.999:
        return 1000.0
    lam = rho * mu

    n_arrivals = int(2.0 * lam * sim_time) + 1000
    interarrivals = np.random.exponential(1 / lam, size=n_arrivals)
    arrival_times = np.cumsum(interarrivals)
    arrival_times = arrival_times[arrival_times <= sim_time]

    if len(arrival_times) == 0:
        return 0.0

    service_times = np.random.exponential(1 / mu, size=len(arrival_times))
    departure_times = np.zeros_like(arrival_times)

    departure_times[0] = arrival_times[0] + service_times[0]
    for i in range(1, len(arrival_times)):
        departure_times[i] = max(arrival_times[i], departure_times[i - 1]) + service_times[i]

    mask = departure_times <= sim_time
    if not np.any(mask):
        return 0.0

    time_in_system = departure_times[mask] - arrival_times[mask]
    return float(np.mean(time_in_system))


def simulate_mm1_algorithm31(rho, mu=1.0, N_min=1000, N_max=200000):
    if rho >= 0.999:
        return 1000.0
    lam = rho * mu

    result = simulate_gi_g_k(
        F_A=lambda: random.expovariate(lam),
        F_B=lambda: random.expovariate(mu),
        lam=lam,
        mean_b=1 / mu,
        K=1,
        N_min=N_min,
        N_max=N_max,
        T_max=0
    )
    return result["E[T]"]


def simulate_md1_algorithm31(rho, mu=1.0, N_min=1000, N_max=200000):
    if rho >= 0.999:
        return 1000.0
    lam = rho * mu
    service_time = 1 / mu

    result = simulate_gi_g_k(
        F_A=lambda: random.expovariate(lam),
        F_B=lambda: service_time,
        lam=lam,
        mean_b=service_time,
        K=1,
        N_min=N_min,
        N_max=N_max,
        T_max=0
    )
    return result["E[T]"]


print("Сравнение методов симуляции M/M/1")

results = []
for i, rho in enumerate(rho_values):
    if (i + 1) % 20 == 0:
        print()

    T_sim_old = simulate_mm1_fast(rho, mu, simulation_time)
    T_sim_alg31 = simulate_mm1_algorithm31(rho, mu, N_min=1000, N_max=100000)
    results.append([rho, T_sim_old, T_sim_alg31])

df = pd.DataFrame(results, columns=["rho", "T_sim_old", "T_sim_alg31"])

df["T_theory"] = 1.0 / (mu * (1 - df["rho"]))
df["E_N_theory"] = df["rho"] / (1 - df["rho"])
df["E_Q_theory"] = df["rho"] ** 2 / (1 - df["rho"])
df["E_W_theory"] = df["rho"] / (mu * (1 - df["rho"]))

rmse_old = np.sqrt(mean_squared_error(df["T_theory"], df["T_sim_old"]))
mae_old = mean_absolute_error(df["T_theory"], df["T_sim_old"])

rmse_alg31 = np.sqrt(mean_squared_error(df["T_theory"], df["T_sim_alg31"]))
mae_alg31 = mean_absolute_error(df["T_theory"], df["T_sim_alg31"])

print(f"\nСтарый метод: RMSE = {rmse_old:.4f}, MAE = {mae_old:.4f}")
print(f"Алгоритм 3.1: RMSE = {rmse_alg31:.4f}, MAE = {mae_alg31:.4f}")

plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(df["rho"], df["T_sim_old"], "b.", alpha=0.4, markersize=4, label="Старый")
plt.plot(df["rho"], df["T_sim_alg31"], "g.", alpha=0.4, markersize=4, label="Алгоритм 3.1")
plt.plot(df["rho"], df["T_theory"], "r-", linewidth=2, label="Теория")
plt.xlabel("rho")
plt.ylabel("E[T]")
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
diff_old = np.abs(df["T_sim_old"] - df["T_theory"])
diff_alg31 = np.abs(df["T_sim_alg31"] - df["T_theory"])
plt.plot(df["rho"], diff_old, "b.", alpha=0.4, markersize=4, label="Старый")
plt.plot(df["rho"], diff_alg31, "g.", alpha=0.4, markersize=4, label="Алгоритм 3.1")
plt.xlabel("rho")
plt.ylabel("Ошибка")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "mm1_methods_comparison.png"), dpi=150, bbox_inches="tight")
plt.show()


print("Сравнение M/M/1 и M/D/1")

results_md1 = []
for i, rho in enumerate(rho_values):
    if (i + 1) % 20 == 0:
        print()

    T_md1 = simulate_md1_algorithm31(rho, mu, N_min=1000, N_max=100000)
    results_md1.append([rho, T_md1])

df_md1 = pd.DataFrame(results_md1, columns=["rho", "T_md1"])

df_md1["T_theory_md1"] = (2 - df_md1["rho"]) / (2 * mu * (1 - df_md1["rho"]))
df_md1["E_N_theory_md1"] = df_md1["rho"] + (df_md1["rho"] ** 2) / (2 * (1 - df_md1["rho"]))
df_md1["E_W_theory_md1"] = df_md1["rho"] / (2 * mu * (1 - df_md1["rho"]))

df_mm1 = df[["rho", "T_sim_alg31", "T_theory"]].copy()
df_comparison = pd.merge(df_mm1, df_md1, on="rho")
df_comparison.rename(columns={"T_theory": "T_theory_mm1"}, inplace=True)

rmse_mm1 = np.sqrt(mean_squared_error(df_comparison["T_theory_mm1"], df_comparison["T_sim_alg31"]))
mae_mm1 = mean_absolute_error(df_comparison["T_theory_mm1"], df_comparison["T_sim_alg31"])
rmse_md1 = np.sqrt(mean_squared_error(df_comparison["T_theory_md1"], df_comparison["T_md1"]))
mae_md1 = mean_absolute_error(df_comparison["T_theory_md1"], df_comparison["T_md1"])

print(f"\nM/M/1: RMSE = {rmse_mm1:.4f}, MAE = {mae_mm1:.4f}")
print(f"M/D/1: RMSE = {rmse_md1:.4f}, MAE = {mae_md1:.4f}")

plt.figure(figsize=(14, 5))
plt.subplot(1, 3, 1)
plt.scatter(df_comparison["rho"], df_comparison["T_sim_alg31"],
            c='blue', alpha=0.5, s=20, label="M/M/1 sim")
plt.plot(df_comparison["rho"], df_comparison["T_theory_mm1"],
         'b-', linewidth=2, label="M/M/1 theory")
plt.scatter(df_comparison["rho"], df_comparison["T_md1"],
            c='orange', alpha=0.5, s=20, label="M/D/1 sim")
plt.plot(df_comparison["rho"], df_comparison["T_theory_md1"],
         'r-', linewidth=2, label="M/D/1 theory")
plt.xlabel("rho")
plt.ylabel("E[T]")
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(1, 3, 2)
high_load = df_comparison[df_comparison["rho"] > 0.7]
if len(high_load) > 0:
    plt.scatter(high_load["rho"], high_load["T_sim_alg31"],
                c='blue', alpha=0.5, s=30, label="M/M/1")
    plt.plot(high_load["rho"], high_load["T_theory_mm1"],
             'b-', linewidth=2)
    plt.scatter(high_load["rho"], high_load["T_md1"],
                c='orange', alpha=0.5, s=30, label="M/D/1")
    plt.plot(high_load["rho"], high_load["T_theory_md1"],
             'r-', linewidth=2)
plt.xlabel("rho")
plt.ylabel("E[T]")
plt.title("rho > 0.7")
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(1, 3, 3)
plt.plot(df_comparison["rho"],
         df_comparison["T_theory_mm1"] - df_comparison["T_theory_md1"],
         'g-', linewidth=2, label="Theory")
plt.plot(df_comparison["rho"],
         df_comparison["T_sim_alg31"] - df_comparison["T_md1"],
         'r--', linewidth=1.5, alpha=0.7, label="Simulation")
plt.xlabel("rho")
plt.ylabel("E[T]_MM1 - E[T]_MD1")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "mm1_md1_comparison.png"), dpi=150, bbox_inches="tight")
plt.show()


print("Анализ точности")

df_comparison["error_mm1"] = np.abs(df_comparison["T_sim_alg31"] - df_comparison["T_theory_mm1"]) / df_comparison["T_theory_mm1"] * 100
df_comparison["error_md1"] = np.abs(df_comparison["T_md1"] - df_comparison["T_theory_md1"]) / df_comparison["T_theory_md1"] * 100

rho_ranges = [(0, 0.3), (0.3, 0.6), (0.6, 0.8), (0.8, 0.95)]
for low, high in rho_ranges:
    mask = (df_comparison["rho"] >= low) & (df_comparison["rho"] < high)
    subset = df_comparison[mask]
    if len(subset) > 0:
        print(f"\nrho [{low:.1f}, {high:.1f}), n={len(subset)}:")
        print(f"  M/M/1 error: {subset['error_mm1'].mean():.2f}%")
        print(f"  M/D/1 error: {subset['error_md1'].mean():.2f}%")

plt.figure(figsize=(10, 5))
plt.plot(df_comparison["rho"], df_comparison["error_mm1"], 'b-', alpha=0.7, label="M/M/1")
plt.plot(df_comparison["rho"], df_comparison["error_md1"], 'r-', alpha=0.7, label="M/D/1")
plt.xlabel("rho")
plt.ylabel("Error %")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(os.path.join(FIG_DIR, "accuracy_analysis.png"), dpi=150, bbox_inches="tight")
plt.show()


print("Отношение M/D/1 к M/M/1")

df_comparison["ratio_theory"] = df_comparison["T_theory_md1"] / df_comparison["T_theory_mm1"]
df_comparison["ratio_sim"] = df_comparison["T_md1"] / df_comparison["T_sim_alg31"]
'''
plt.figure(figsize=(10, 5))
plt.plot(df_comparison["rho"], df_comparison["ratio_theory"], 'b-', linewidth=2, label="Theory")
plt.scatter(df_comparison["rho"], df_comparison["ratio_sim"], c='red', alpha=0.5, s=20, label="Simulation")
plt.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
plt.xlabel("rho")
plt.ylabel("E[T]_MD1 / E[T]_MM1")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
'''


print("Теоретические формулы")

print("\nM/M/1:")
print("  E[N] = rho/(1-rho)")
print("  E[T] = 1/(1-rho)")
print("  E[W] = rho/(1-rho)")
print("  E[Q] = rho^2/(1-rho)")

print("\nM/D/1:")
print("  E[N] = rho + rho^2/(2(1-rho))")
print("  E[T] = (2-rho)/(2(1-rho))")
print("  E[W] = rho/(2(1-rho))")
print("  E[Q] = rho^2/(2(1-rho))")


print("Итоги")

print(f"Точек: {len(df)}")
print(f"Диапазон rho: [{rho_values[0]:.3f}, {rho_values[-1]:.3f}]")
print(f"\nСредняя ошибка M/M/1: {df_comparison['error_mm1'].mean():.2f}%")
print(f"Средняя ошибка M/D/1: {df_comparison['error_md1'].mean():.2f}%")

idx_05 = (np.abs(df_comparison["rho"] - 0.5)).argmin()
idx_09 = (np.abs(df_comparison["rho"] - 0.9)).argmin()

print(f"\nВыигрыш M/D/1:")
print(f"  rho={df_comparison.iloc[idx_05]['rho']:.2f}: {100 * (1 - df_comparison.iloc[idx_05]['ratio_theory']):.1f}%")
print(f"  rho={df_comparison.iloc[idx_09]['rho']:.2f}: {100 * (1 - df_comparison.iloc[idx_09]['ratio_theory']):.1f}%")
