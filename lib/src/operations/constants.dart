import 'dart:math';

/// 2π
const double twoPi = 2.0 * pi;

/// π/2
const double halfPi = 0.5 * pi;

/// Convert degrees to radians.
const double deg2rad = pi / 180.0;

/// Convert radians to degrees.
const double rad2deg = 180.0 / pi;

/// Convert seconds to degrees.
const double sec2deg = 1.0 / 60.0 / 60.0;

/// Convert seconds to days.
const double sec2day = sec2deg / 24.0;

/// Convert arcseconds to radians.
const double asec2rad = sec2deg * deg2rad;

/// Convert ten-thousandths of an arcsecond to radians.
const double ttasec2rad = asec2rad / 10000.0;

/// Convert milliarcseconds to radians.
const double masec2rad = asec2rad / 1000.0;

/// Astronomical unit _(km)_.
const double astronomicalUnit = 149597870.0;

/// Convert milliseconds to seconds.
const double msec2sec = 1e-3;

/// Speed of light _(km/s)_.
const double speedOfLight = 299792.458;

/// Seconds per day.
const double secondsPerDay = 86400.0;

/// Convert seconds to minutes.
const double sec2min = 1.0 / 60.0;

/// Seconds per sidereal day.
const double secondsPerSiderealDay = 86164.0905;

/// Seconds per week.
const double secondsPerWeek = secondsPerDay * 7.0;
