"""D2Q9-TRT chi_kappa-to-ce thermal-model regression checks.

These compact tests exercise the algebra used by the Fortran implementation and
two thermal-only limits: periodic sinusoidal diffusion and halfway-ABB steady
conduction.  They intentionally do not exercise the coupled OpenACC RB solver.
"""

from __future__ import annotations

import math

import numpy as np


CS2 = 1.0 / 3.0
EX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=int)
EY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=int)
OPPOSITE = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=int)
OMEGA = np.array([4.0 / 9.0] + [1.0 / 9.0] * 4 + [1.0 / 36.0] * 4)
LAMBDA = np.array([-5.0 / 9.0] + [1.0 / 9.0] * 4 + [1.0 / 36.0] * 4)
QK = 3.0 - math.sqrt(3.0)
QE = 4.0 * math.sqrt(3.0) - 6.0
QNU = 4.0 * math.sqrt(3.0) - 6.0

M = np.array(
    [
        [1, 1, 1, 1, 1, 1, 1, 1, 1],
        [-4, -1, -1, -1, -1, 2, 2, 2, 2],
        [4, -2, -2, -2, -2, 1, 1, 1, 1],
        [0, 1, 0, -1, 0, 1, -1, -1, 1],
        [0, -2, 0, 2, 0, 1, -1, -1, 1],
        [0, 0, 1, 0, -1, 1, 1, -1, -1],
        [0, 0, -2, 0, 2, 1, 1, -1, -1],
        [0, 1, -1, 1, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, -1, 1, -1],
    ],
    dtype=float,
)


def ce_from_para_a(para_a: float) -> float:
    return (para_a + 4.0) / 6.0


def ce_from_chi_kappa(chi_kappa: float, diffusivity: float) -> float:
    return diffusivity / ((1.0 - chi_kappa) * (1.0 / QK - 0.5))


def para_a_from_chi_kappa(chi_kappa: float, diffusivity: float) -> float:
    return 6.0 * ce_from_chi_kappa(chi_kappa, diffusivity) - 4.0


def chi_kappa_is_valid(chi_kappa: float, diffusivity: float) -> bool:
    chi_kappa_max = 1.0 - diffusivity / ((3.0 / 5.0) * (1.0 / QK - 0.5))
    return 0.0 <= chi_kappa < chi_kappa_max


def thermal_weights(para_a: float) -> np.ndarray:
    return np.array(
        [-(5.0 * para_a + 2.0) / 18.0]
        + [(para_a + 4.0) / 18.0] * 4
        + [(para_a + 4.0) / 72.0] * 4
    )


def equilibrium(
    temperature: np.ndarray,
    para_a: float,
    velocity_x: float | np.ndarray = 0.0,
    velocity_y: float | np.ndarray = 0.0,
    pressure: float | np.ndarray = 0.0,
) -> np.ndarray:
    ce = ce_from_para_a(para_a)
    weights = thermal_weights(para_a)
    speed_squared = velocity_x * velocity_x + velocity_y * velocity_y
    result = np.empty((9,) + temperature.shape, dtype=float)
    for alpha in range(9):
        eu = EX[alpha] * velocity_x + EY[alpha] * velocity_y
        result[alpha] = (
            weights[alpha] * temperature
            + OMEGA[alpha]
            * temperature
            * (eu / CS2 + 0.5 * eu * eu / (CS2 * CS2) - 0.5 * speed_squared / CS2)
            + LAMBDA[alpha] * temperature * pressure / CS2
        )
    return result


def collide_pure_diffusion(g: np.ndarray, chi_kappa: float, diffusivity: float) -> np.ndarray:
    ce = ce_from_chi_kappa(chi_kappa, diffusivity)
    para_a = 6.0 * ce - 4.0
    temperature = np.sum(g, axis=0)
    geq = equilibrium(temperature, para_a)
    neq = g - geq
    gradient_x = -2.0 * np.sum(EX[:, None, None] * neq, axis=0) / (2.0 * diffusivity + ce)
    gradient_y = -2.0 * np.sum(EY[:, None, None] * neq, axis=0) / (2.0 * diffusivity + ce)
    source_x = chi_kappa * ce * gradient_x
    source_y = chi_kappa * ce * gradient_y
    moments = np.einsum("ij,j...->i...", M, g)
    equilibrium_moments = np.einsum("ij,j...->i...", M, geq)
    source_moments = np.zeros_like(moments)
    source_moments[3] = source_x
    source_moments[4] = -source_x
    source_moments[5] = source_y
    source_moments[6] = -source_y
    rates = np.array([0.0, QE, QE, QK, QK, QK, QK, QNU, QNU])[:, None, None]
    post_moments = (
        moments
        - rates * (moments - equilibrium_moments)
        + (1.0 - 0.5 * rates) * source_moments
    )
    return np.einsum("ij,j...->i...", np.linalg.inv(M), post_moments)


def stream_periodic(post: np.ndarray) -> np.ndarray:
    streamed = np.empty_like(post)
    for alpha in range(9):
        streamed[alpha] = np.roll(post[alpha], shift=(EY[alpha], EX[alpha]), axis=(0, 1))
    return streamed


def test_algebra() -> None:
    para_a = -1.54
    ce = ce_from_para_a(para_a)
    temperature = np.array([[0.73]])
    velocity_x = 0.031
    velocity_y = -0.027
    pressure = 2.5e-4
    geq = equilibrium(temperature, para_a, velocity_x, velocity_y, pressure)[:, 0, 0]
    moments = M @ geq
    speed_squared = velocity_x**2 + velocity_y**2
    expected = np.array(
        [
            temperature.item(),
            (para_a + 3.0 * speed_squared + 6.0 * pressure) * temperature.item(),
            (-(3.0 * para_a + 4.0) / 2.0 - 3.0 * speed_squared - 9.0 * pressure)
            * temperature.item(),
            velocity_x * temperature.item(),
            -velocity_x * temperature.item(),
            velocity_y * temperature.item(),
            -velocity_y * temperature.item(),
            (velocity_x**2 - velocity_y**2) * temperature.item(),
            velocity_x * velocity_y * temperature.item(),
        ]
    )
    assert np.max(np.abs(moments - expected)) < 1.0e-13

    source_vector = np.array([0.017, -0.023])
    source = OMEGA * (EX * source_vector[0] + EY * source_vector[1]) / CS2
    expected_source_moments = np.array(
        [0.0, 0.0, 0.0, source_vector[0], -source_vector[0], source_vector[1], -source_vector[1], 0.0, 0.0]
    )
    assert np.max(np.abs(M @ source - expected_source_moments)) < 1.0e-14

    weights = thermal_weights(para_a)
    assert abs(np.sum(weights) - 1.0) < 1.0e-15
    assert abs(2.0 * weights[1] + 4.0 * weights[5] - ce) < 1.0e-15
    standard_para_a = -2.0
    standard_ce = ce_from_para_a(standard_para_a)
    assert abs(standard_ce - CS2) < 1.0e-15
    assert np.max(np.abs(thermal_weights(standard_para_a) - OMEGA)) < 1.0e-15
    assert abs((6.0 * standard_ce - 4.0) - (-2.0)) < 1.0e-15
    assert abs((4.0 - 9.0 * standard_ce) - 1.0) < 1.0e-15
    assert abs((1.0 / QK - 0.5) - 1.0 / math.sqrt(12.0)) < 1.0e-15
    assert abs((1.0 / QE - 0.5) - 1.0 / math.sqrt(3.0)) < 1.0e-15
    assert abs(QE - QNU) < 1.0e-15

    rng = np.random.default_rng(20260818)
    g = rng.random(9)
    geq = rng.random(9)
    geq *= np.sum(g) / np.sum(geq)
    source = OMEGA * (EX * source_vector[0] + EY * source_vector[1]) / CS2
    direct = (
        g
        - QE * (0.5 * (g + g[OPPOSITE]) - 0.5 * (geq + geq[OPPOSITE]))
        - QK * (0.5 * (g - g[OPPOSITE]) - 0.5 * (geq - geq[OPPOSITE]))
        + (1.0 - 0.5 * QK) * source
    )
    relaxation = np.diag([0.0, QE, QE, QK, QK, QK, QK, QNU, QNU])
    moment_space = np.linalg.solve(
        M,
        M @ g - relaxation @ (M @ g - M @ geq) + (np.eye(9) - 0.5 * relaxation) @ (M @ source),
    )
    assert np.max(np.abs(direct - moment_space)) < 2.0e-14


def test_parameter_bounds() -> None:
    diffusivity = 0.02
    chi_kappa = 0.5
    ce = ce_from_chi_kappa(chi_kappa, diffusivity)
    para_a = para_a_from_chi_kappa(chi_kappa, diffusivity)
    reconstructed = (1.0 - chi_kappa) * ce * (1.0 / QK - 0.5)
    assert chi_kappa_is_valid(chi_kappa, diffusivity)
    assert abs(para_a - (6.0 * ce - 4.0)) < 1.0e-15
    assert abs(reconstructed - diffusivity) < 1.0e-15

    chi_kappa_max = 1.0 - diffusivity / ((3.0 / 5.0) * (1.0 / QK - 0.5))
    assert chi_kappa_is_valid(0.0, diffusivity)
    assert not chi_kappa_is_valid(np.nextafter(0.0, -math.inf), diffusivity)
    assert not chi_kappa_is_valid(chi_kappa_max, diffusivity)
    assert not chi_kappa_is_valid(1.0, diffusivity)


def test_periodic_diffusion() -> list[tuple[float, float, float, float, float]]:
    nx, ny = 128, 4
    diffusivity = 0.02
    wave_number = 2.0 * math.pi / nx
    x = np.arange(nx, dtype=float)
    mode = np.cos(wave_number * x)[None, :]
    initial_temperature = np.ones((ny, nx)) + 1.0e-3 * mode
    results: list[tuple[float, float, float, float, float]] = []
    for chi_kappa in (0.0, 0.5, 0.8):
        ce = ce_from_chi_kappa(chi_kappa, diffusivity)
        para_a = 6.0 * ce - 4.0
        g = equilibrium(initial_temperature, para_a)
        amplitudes: dict[int, float] = {}
        for step in range(1601):
            if step in (400, 1600):
                temperature = np.sum(g, axis=0)
                amplitudes[step] = 2.0 * np.mean((temperature - np.mean(temperature)) * mode)
            if step < 1600:
                g = stream_periodic(collide_pure_diffusion(g, chi_kappa, diffusivity))
        measured = -math.log(amplitudes[1600] / amplitudes[400]) / (wave_number**2 * 1200.0)
        relative_error = abs(measured - diffusivity) / diffusivity
        assert relative_error < 5.0e-3
        results.append((chi_kappa, para_a, ce, measured, relative_error))
    return results


def test_halfway_abb_conduction() -> list[tuple[float, float, float, float]]:
    nx, ny = 8, 48
    hot, cold = 0.5, -0.5
    diffusivity = 0.02
    y = (np.arange(ny, dtype=float) + 0.5) / ny
    exact = hot + y[:, None] * (cold - hot) + np.zeros((ny, nx))
    results: list[tuple[float, float, float, float]] = []
    for chi_kappa in (0.0, 0.5, 0.8):
        ce = ce_from_chi_kappa(chi_kappa, diffusivity)
        para_a = 6.0 * ce - 4.0
        weights = thermal_weights(para_a)
        g = equilibrium(exact, para_a)
        for _ in range(200):
            post = collide_pure_diffusion(g, chi_kappa, diffusivity)
            g = stream_periodic(post)
            for alpha in (2, 5, 6):
                g[alpha, 0, :] = -post[OPPOSITE[alpha], 0, :] + 2.0 * weights[alpha] * hot
            for alpha in (4, 7, 8):
                g[alpha, -1, :] = -post[OPPOSITE[alpha], -1, :] + 2.0 * weights[alpha] * cold
        maximum_error = float(np.max(np.abs(np.sum(g, axis=0) - exact)))
        assert maximum_error < 1.0e-12
        results.append((chi_kappa, para_a, ce, maximum_error))
    return results


def main() -> None:
    test_algebra()
    test_parameter_bounds()
    diffusion_results = test_periodic_diffusion()
    conduction_results = test_halfway_abb_conduction()
    print("algebra: PASS")
    print("chi_kappa bounds and chi_kappa-to-ce kappa reconstruction: PASS")
    for chi_kappa, para_a, ce, measured, relative_error in diffusion_results:
        print(
            f"periodic diffusion chi_kappa={chi_kappa:.12g}, paraA={para_a:.12g}, ce={ce:.12g}: "
            f"kappa={measured:.12g}, relative_error={relative_error:.3e}"
        )
    for chi_kappa, para_a, ce, maximum_error in conduction_results:
        print(
            f"halfway ABB conduction chi_kappa={chi_kappa:.12g}, paraA={para_a:.12g}, ce={ce:.12g}: "
            f"max_error={maximum_error:.3e}"
        )


if __name__ == "__main__":
    main()
