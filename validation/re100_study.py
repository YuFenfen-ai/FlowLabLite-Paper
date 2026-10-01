"""Compatibility wrapper for the generalized Ghia validation study."""

from validation.ghia_study import *  # noqa: F401,F403
from validation.ghia_study import _is_time_converged, main


if __name__ == "__main__":
    raise SystemExit(main())