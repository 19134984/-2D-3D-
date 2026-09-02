import unittest


def restart_time_tolerance(time_unit, restart_time, sample_interval):
    scale_tolerance = 1.0e-10 * max(1.0, abs(restart_time), abs(sample_interval))
    return max(scale_tolerance, 0.51 / time_unit)


class RestartTimeToleranceTest(unittest.TestCase):
    def test_integer_lattice_clock_roundoff_is_accepted(self):
        time_unit = 35472.40053901060
        scheduled_time = 100.0
        restart_itc = round(scheduled_time * time_unit)
        restart_time = restart_itc / time_unit

        old_tolerance = 1.0e-10 * max(1.0, abs(restart_time), 1.0)
        mismatch = abs(scheduled_time - restart_time)

        self.assertGreater(mismatch, old_tolerance)
        self.assertLessEqual(
            mismatch,
            restart_time_tolerance(time_unit, restart_time, 1.0),
        )

    def test_one_full_lattice_step_is_still_rejected(self):
        time_unit = 35472.40053901060
        tolerance = restart_time_tolerance(time_unit, 100.0, 1.0)

        self.assertGreater(1.0 / time_unit, tolerance)


if __name__ == "__main__":
    unittest.main()
