"""Adaptive broadening and wavelength-shift mixin for forward models."""

import typing as t

import numpy as np
import numpy.typing as npt

from taurex.binning import FluxBinner
from taurex.core import fitparam
from taurex.util import create_grid_res

from .core import ForwardModelMixin


class AdaptiveForwardModel(ForwardModelMixin):
    """Adds wavelength shifting and low-resolution broadening to a model."""

    def __init_mixin__(
        self,
        wlshift: t.Optional[t.Sequence[float]] = None,
        wlbroadening_method: str = "box",
        wlbroadening_width: float = 3,
        wlbroadening_width2: float = 0,
        wlbroadening_default: float = 0,
        max_wlbroadening: float = 0.1,
        factor_cut: int = 5,
        wlres: float = 15000,
    ) -> None:
        self._wlshift = list(wlshift or [1.0, 0.0])
        self._wlbroadening_method = wlbroadening_method
        self._wlbroadening_width = wlbroadening_width
        self._wlbroadening_width2 = wlbroadening_width2
        self._wlbroadening_default = wlbroadening_default
        self._max_wlbroadening = max_wlbroadening
        self._factor_cut = factor_cut
        self._wlres = wlres
        self.generate_waveshift_fitting_params()

    @staticmethod
    def gaussian(
        x: npt.NDArray[np.float64], mean: float, std: float
    ) -> npt.NDArray[np.float64]:
        return (
            1.0
            / (np.sqrt(2.0 * np.pi) * std)
            * np.exp(-np.power((x - mean) / std, 2.0) / 2.0)
        )

    def box_moving_average(
        self,
        output: t.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], t.Any, t.Any],
    ) -> t.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        windows = max(1, int(self._wlbroadening_width))
        kernel = np.ones(windows) / windows
        return output[0], np.convolve(output[1], kernel, mode="same")

    def gaussian_moving_average(
        self,
        output: t.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], t.Any, t.Any],
    ) -> t.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        windows = max(1, int(self._wlbroadening_width))
        x = np.linspace(-5 * windows, 5 * windows, windows)
        kernel = self.gaussian(x, 0.0, windows)
        kernel = kernel / np.sum(kernel)
        return output[0], np.convolve(output[1], kernel, mode="same")

    def low_res_convolved(
        self,
        output: t.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], t.Any, t.Any],
    ) -> t.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        native_grid, final_flux, _, _ = output
        new_grid = create_grid_res(self._wlres, native_grid[0], native_grid[-1])
        flux_binner = FluxBinner(10000.0 / new_grid[:, 0], new_grid[:, 1])
        binned = flux_binner.bindown(native_grid, final_flux)
        wavelengths = 10000.0 / binned[0]

        convolved = np.zeros(len(wavelengths))
        for index, wavelength in enumerate(wavelengths):
            std = np.clip(
                self._wlbroadening_default
                + self._wlbroadening_width * wavelength
                + self._wlbroadening_width2 * wavelength**2,
                a_min=1e-20,
                a_max=self._max_wlbroadening,
            )
            if index != len(wavelengths) - 1:
                spacing = np.abs(wavelengths[index + 1] - wavelength)
            else:
                spacing = np.abs(wavelength - wavelengths[index - 1])
            window = max(1, int(self._factor_cut * std / spacing))
            start = max(0, index - window)
            stop = min(len(wavelengths), index + window + 1)
            x = wavelengths[start:stop]
            y = self.gaussian(x, wavelength, std)
            convolved[start:stop] += binned[1][index] * y / np.sum(y)

        return 10000.0 / wavelengths, convolved

    def low_res_convolved_r(
        self,
        output: t.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], t.Any, t.Any],
    ) -> t.Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
        native_grid, final_flux, _, _ = output
        new_grid = create_grid_res(self._wlres, native_grid[0], native_grid[-1])
        flux_binner = FluxBinner(10000.0 / new_grid[:, 0], new_grid[:, 1])
        binned = flux_binner.bindown(native_grid, final_flux)
        wavelengths = 10000.0 / binned[0]

        convolved = np.zeros(len(wavelengths))
        for index, wavelength in enumerate(wavelengths):
            std = np.clip(
                0.5
                * wavelength
                / (
                    self._wlbroadening_default
                    + self._wlbroadening_width * wavelength
                    + self._wlbroadening_width2 * wavelength**2
                ),
                a_min=1e-20,
                a_max=self._max_wlbroadening,
            )
            if index != len(wavelengths) - 1:
                spacing = np.abs(wavelengths[index + 1] - wavelength)
            else:
                spacing = np.abs(wavelength - wavelengths[index - 1])
            window = max(1, int(self._factor_cut * std / spacing))
            start = max(0, index - window)
            stop = min(len(wavelengths), index + window + 1)
            x = wavelengths[start:stop]
            y = self.gaussian(x, wavelength, std)
            convolved[start:stop] += binned[1][index] * y / np.sum(y)

        return 10000.0 / wavelengths, convolved

    def model(
        self,
        wngrid: t.Optional[npt.NDArray[np.float64]] = None,
        cutoff_grid: t.Optional[bool] = True,
    ) -> t.Tuple[
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
        t.Optional[npt.NDArray[np.float64]],
        t.Optional[t.Any],
    ]:
        output = super().model(wngrid=wngrid, cutoff_grid=cutoff_grid)

        if self._wlbroadening_method == "binned_convolution":
            adapted_grid, adapted_flux = self.low_res_convolved(output)
        elif self._wlbroadening_method == "binned_convolution_R":
            adapted_grid, adapted_flux = self.low_res_convolved_r(output)
        elif self._wlbroadening_method == "gaussian_box" and self._wlbroadening_width > 0:
            adapted_grid, adapted_flux = self.gaussian_moving_average(output)
        elif self._wlbroadening_method == "box" and self._wlbroadening_width > 0:
            adapted_grid, adapted_flux = self.box_moving_average(output)
        else:
            adapted_grid = output[0]
            adapted_flux = output[1]

        wavelengths = 10000.0 / adapted_grid[::-1]
        shifted_grid = 10000.0 / np.polyval(self._wlshift, wavelengths)[::-1]
        return shifted_grid, adapted_flux, output[2], None

    def generate_waveshift_fitting_params(self) -> None:
        bounds = (-1.0, 1.0)
        for index in range(len(self._wlshift)):
            num = index + 1

            def read_shift(self, index=index):
                return self._wlshift[index]

            def write_shift(self, value, index=index):
                self._wlshift[index] = value

            self.add_fittable_param(
                f"wlshift_{num}",
                f"$wlshift_{num}$",
                read_shift,
                write_shift,
                "linear",
                False,
                bounds,
            )

    @fitparam(
        param_name="Wbroad",
        param_latex="$Wbroad$",
        default_fit=False,
        default_bounds=[0.0, 0.0002],
    )
    def waveBroad(self) -> float:  # noqa: N802
        return self._wlbroadening_width

    @waveBroad.setter
    def waveBroad(self, value: float) -> None:  # noqa: N802
        self._wlbroadening_width = value

    @fitparam(
        param_name="Wbroad2",
        param_latex="$Wbroad2$",
        default_fit=False,
        default_bounds=[0.0, 0.0002],
    )
    def waveBroad2(self) -> float:  # noqa: N802
        return self._wlbroadening_width2

    @waveBroad2.setter
    def waveBroad2(self, value: float) -> None:  # noqa: N802
        self._wlbroadening_width2 = value

    @fitparam(
        param_name="Dbroad",
        param_latex="$Dbroad$",
        default_fit=False,
        default_bounds=[0.0, 0.0002],
    )
    def defaultWaveBroad(self) -> float:  # noqa: N802
        return self._wlbroadening_default

    @defaultWaveBroad.setter
    def defaultWaveBroad(self, value: float) -> None:  # noqa: N802
        self._wlbroadening_default = value

    def write(self, output):
        model = super().write(output)
        model.write_scalar("Wbroad", self._wlbroadening_width)
        model.write_scalar("Dbroad", self._wlbroadening_default)
        model.write_array("Wshift", np.array(self._wlshift))
        model.write_string("Wbroad_method", self._wlbroadening_method)
        return model

    @classmethod
    def input_keywords(cls) -> t.Tuple[str, ...]:
        return ("adaptivemodel", "adaptive")