"""Сравнение среднего числа заявок E[N] в системах M/M/1 и M/D/1.

Используются три способа расчёта:
  1. "Быстрая" симуляция M/M/1 на numpy (векторизованное построение
     потока поступлений и процесса обслуживания);
  2. Алгоритм 3.1 моделирования системы GI/G/K из книги
     П. П. Бочарова, А. В. Печинкина «Теория массового обслуживания»,
     стр. 110–114;
  3. Аналитические формулы для E[N] в системах M/M/1 и M/D/1 (стр. 114).

Для каждой загрузки rho из [0.01, 0.95] выводится таблица значений E[N],
в конце — среднее время работы методов и RMSE относительно теории.
График сохраняется в figures/mm1_md1_mean_number.png.
"""

import os
import time
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error

FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

np.random.seed(42)
random.seed(42)

mu = 1.0
num_rho = 50
rho_values = np.linspace(0.01, 0.95, num_rho)
simulation_time = 500
N_min = 1000
N_max = 200000


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
        x_a -= T_x
        x_b = [xb - T_x for xb in x_b]

        if N_a >= N_min:
            T_n += n * T_x
            T_s += T_x

        if x_a <= eps:
            N_a += 1
            n += 1
            x_a = F_A()
            if n <= K:
                x_b.append(F_B())
        elif min(x_b) <= eps:
            j = x_b.index(min(x_b))
            x_b.pop(j)
            n -= 1
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
    return {"E[N]": E_N, "E[T]": E_T, "E[W]": E_W, "E[Q]": E_Q}


# "Быстрая" симуляция M/M/1 (построение траектории на numpy)
def simulate_mm1_fast(rho, mu=1.0, sim_time=500):
    if rho >= 0.999:
        return float("inf")
    lam = rho * mu
    n_arrivals = int(2 * lam * sim_time) + 1000
    interarrivals = np.random.exponential(1 / lam, size=n_arrivals)
    arrival_times = np.cumsum(interarrivals)
    arrival_times = arrival_times[arrival_times <= sim_time]
    if len(arrival_times) == 0:
        return 0.0
    service_times = np.random.exponential(1 / mu, size=len(arrival_times))
    depart = np.zeros_like(arrival_times)
    depart[0] = arrival_times[0] + service_times[0]
    for i in range(1, len(arrival_times)):
        depart[i] = max(arrival_times[i], depart[i - 1]) + service_times[i]
    mask = depart <= sim_time
    if not np.any(mask):
        return 0.0
    time_in_sys = depart[mask] - arrival_times[mask]
    N = lam * np.mean(time_in_sys)
    return N


def simulate_mm1_alg31(rho, mu=1.0, N_min=1000, N_max=200000):
    if rho >= 0.999:
        return float("inf")
    lam = rho * mu
    res = simulate_gi_g_k(
        F_A=lambda: random.expovariate(lam),
        F_B=lambda: random.expovariate(mu),
        lam=lam, mean_b=1 / mu, K=1,
        N_min=N_min, N_max=N_max, T_max=0,
    )
    return res["E[N]"]


def simulate_md1_alg31(rho, mu=1.0, N_min=1000, N_max=200000):
    if rho >= 0.999:
        return float("inf")
    lam = rho * mu
    service_time = 1 / mu
    res = simulate_gi_g_k(
        F_A=lambda: random.expovariate(lam),
        F_B=lambda: service_time,
        lam=lam, mean_b=service_time, K=1,
        N_min=N_min, N_max=N_max, T_max=0,
    )
    return res["E[N]"]


# Аналитические формулы из книги (стр. 114)
def N_mm1_theory(rho):
    return rho / (1 - rho)


def N_md1_theory(rho):
    return rho + (rho ** 2) / (2 * (1 - rho))


# Основной расчёт
results = []
times = []

print("rho\tN_mm1_fast\tN_mm1_alg31\tN_mm1_theory\tN_md1_alg31\tN_md1_theory")
for rho in rho_values:
    # M/M/1 fast
    t0 = time.time()
    n_mm1_fast = simulate_mm1_fast(rho, mu, simulation_time)
    t1 = time.time()
    # M/M/1 alg31
    n_mm1_alg31 = simulate_mm1_alg31(rho, mu, N_min, N_max)
    t2 = time.time()
    # M/D/1 alg31
    n_md1_alg31 = simulate_md1_alg31(rho, mu, N_min, N_max)
    t3 = time.time()
    # theory
    n_mm1_th = N_mm1_theory(rho)
    n_md1_th = N_md1_theory(rho)

    results.append([rho, n_mm1_fast, n_mm1_alg31, n_mm1_th, n_md1_alg31, n_md1_th])
    times.append([t1 - t0, t2 - t1, t3 - t2])
    print(f"{rho:.4f}\t{n_mm1_fast:.3f}\t\t{n_mm1_alg31:.3f}\t\t{n_mm1_th:.3f}\t\t{n_md1_alg31:.3f}\t\t{n_md1_th:.3f}")

# Среднее время работы
df_times = pd.DataFrame(times, columns=["mm1_fast", "mm1_alg31", "md1_alg31"])
print("\nСреднее время работы на один rho (сек):")
print(f"  M/M/1 fast   : {df_times['mm1_fast'].mean():.4f}")
print(f"  M/M/1 alg31  : {df_times['mm1_alg31'].mean():.4f}")
print(f"  M/D/1 alg31  : {df_times['md1_alg31'].mean():.4f}")

# Ошибка RMSE относительно теории
df = pd.DataFrame(results, columns=["rho", "N_mm1_fast", "N_mm1_alg31", "N_mm1_th", "N_md1_alg31", "N_md1_th"])
rmse_fast = np.sqrt(mean_squared_error(df["N_mm1_th"], df["N_mm1_fast"]))
rmse_alg31 = np.sqrt(mean_squared_error(df["N_mm1_th"], df["N_mm1_alg31"]))
rmse_md1 = np.sqrt(mean_squared_error(df["N_md1_th"], df["N_md1_alg31"]))
print("\nRMSE (относительно теории):")
print(f"  M/M/1 fast   : {rmse_fast:.4f}")
print(f"  M/M/1 alg31  : {rmse_alg31:.4f}")
print(f"  M/D/1 alg31  : {rmse_md1:.4f}")

# График
plt.figure(figsize=(10, 5))
plt.plot(df["rho"], df["N_mm1_th"], "k-", label="M/M/1 theory")
plt.plot(df["rho"], df["N_mm1_alg31"], "r.", markersize=2, label="M/M/1 sim")
plt.plot(df["rho"], df["N_md1_th"], "b-", label="M/D/1 theory")
plt.plot(df["rho"], df["N_md1_alg31"], "g.", markersize=2, label="M/D/1 sim")
plt.xlabel("rho")
plt.ylabel("N")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(os.path.join(FIG_DIR, "mm1_md1_mean_number.png"), dpi=150, bbox_inches="tight")
plt.show()
