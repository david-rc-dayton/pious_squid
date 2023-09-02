import 'dart:math';

import 'package:pious_squid/pious_squid.dart';

/// Convert radians to degrees.
const rad2deg = 180 / pi;

void main() {
  // Define some values for future calculations.
  final epoch = EpochUTC.fromDate(2016, 12, 12, 21, 35, 15.425);
  final location = Geodetic.fromDegrees(38.4, -104.2, 1.7);
  final position = location.toITRF(epoch).toJ2000().position;

  // ---- Earth -----------------------------------------------------------------

  // Calculate orbit mean motion (rad/s) from semimajor axis (km).
  final meanMotion = Earth.smaToMeanMotion(42164.1);
  print(meanMotion); // => 0.00007292133916998565 (rad)
  print(meanMotion * rad2deg); // => 0.004178084970882192 (deg)

  // Calculate orbit semimajor axis (km) given revolutions per day.
  final semimajorAxis = Earth.revsPerDayToSma(
    15.501, // revs/day
  );
  print(semimajorAxis); // => 6794.570831048931 (km)

  // ---- Sun -----------------------------------------------------------------

  // Get Sun inertial position (km).
  final sunPosition = Sun.position(epoch);
  print(sunPosition);
  // => [-22866135.59362657, -133484239.13875942, -57866314.76604423]

  // Get Sun apparent inertial position as percieved from Earth (km), since
  // light takes a bit of time to travel.
  final sunApparentPosition = Sun.positionApparent(epoch);
  print(sunApparentPosition);
  // => [-22880829.99855808, -133482204.21388586, -57865432.56112264]

  // Get Sun angular diameter from an inertial observer position (rad).
  final sunDiameter = Sun.diameter(position, sunPosition);
  print(sunDiameter); // => 0.009445185896037878 (rad)
  print(sunDiameter * rad2deg); // => 0.5411692885594611 (deg)

  // Check if the observer is eclipsed by the Earth.
  final shadow = Sun.shadow(epoch, position);
  print(shadow); // => false (not in shadow)

  // Calculate the lighting ratio of the observer (ratio).
  final lightingRatio = Sun.lightingRatio(position, sunPosition);
  print(lightingRatio); // => 1.0 (fully illuminated)

  // ---- Moon ----------------------------------------------------------------

  // Get Moon inertial position (km).
  final moonPosition = Moon.position(epoch);
  print(moonPosition);
  // => [148559.24735268456, 311217.93470951146, 100876.67340196588]

  // Get Moon angular diameter from an observer position (rad).
  final moonDiameter = Moon.diameter(position, moonPosition);
  print(moonDiameter); // => 0.009632157811104816 (rad)
  print(moonDiameter * rad2deg); // => 0.5518819901802752 (deg)

  // Calculate Moon illumination.
  final moonIllumination = Moon.illumination(epoch, position);
  print(moonIllumination); // => 0.980361106792833 (nearly fully illuminated)
}
