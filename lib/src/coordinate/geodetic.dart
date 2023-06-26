import 'dart:math';

import 'package:pious_squid/src/body/body_base.dart';
import 'package:pious_squid/src/coordinate/coordinate_base.dart';
import 'package:pious_squid/src/operations/constants.dart';
import 'package:pious_squid/src/operations/functions.dart';
import 'package:pious_squid/src/operations/operations_base.dart';
import 'package:pious_squid/src/time/time_base.dart';

/// Geodetic coordinates.
class Geodetic {
  /// Create a new [Geodetic] object.
  Geodetic(this.latitude, this.longitude, this.altitude);

  /// Create a new [Geodetic] object using degrees for angular fields.
  factory Geodetic.fromDegrees(
          final double latDeg, final double lonDeg, final double alt) =>
      Geodetic(latDeg * deg2rad, lonDeg * deg2rad, alt);

  /// Latitude _(rad)_.
  final double latitude;

  /// Longitude _(rad)_.
  final double longitude;

  /// Altitude _(km)_.
  final double altitude;

  /// String representation of this object.
  @override
  String toString() => [
        '[Geodetic]',
        '  Latitude:  ${latitudeDegrees.toStringAsFixed(6)}°',
        '  Longitude: ${longitudeDegrees.toStringAsFixed(6)}°',
        '  Altitude:  ${altitude.toStringAsFixed(3)} km',
      ].join('\n');

  /// Latitude _(deg)_.
  double get latitudeDegrees => latitude * rad2deg;

  /// Longitude _(deg)_.
  double get longitudeDegrees => longitude * rad2deg;

  /// Convert this to [ITRF] coordinates at the provided UTC [epoch].
  ITRF toITRF(final EpochUTC epoch) {
    final sLat = sin(latitude);
    final cLat = cos(latitude);
    final nVal =
        Earth.radiusEquator / sqrt(1 - Earth.eccentricitySquared * sLat * sLat);
    final r = Vector3D(
        (nVal + altitude) * cLat * cos(longitude),
        (nVal + altitude) * cLat * sin(longitude),
        (nVal * (1 - Earth.eccentricitySquared) + altitude) * sLat);
    return ITRF(epoch, r, Vector3D.origin);
  }

  /// Calculate the angular distance _(rad)_ between this and another
  /// [Geodetic] object.
  double angle(final Geodetic g,
          {final AngularDistanceMethod method =
              AngularDistanceMethod.haversine}) =>
      angularDistance(longitude, latitude, g.longitude, g.latitude,
          method: method);

  /// Calculate the angular distance _(°)_ between this and another
  /// [Geodetic] object.
  double angleDegrees(final Geodetic g,
          {final AngularDistanceMethod method =
              AngularDistanceMethod.haversine}) =>
      angle(g, method: method) * rad2deg;

  /// Calculate the distance _(km)_ along the Earth's surface between this and
  /// another [Geodetic] object.
  double distance(final Geodetic g,
          {final AngularDistanceMethod method =
              AngularDistanceMethod.haversine}) =>
      angle(g, method: method) * Earth.radiusMean;

  /// Calculate the angular field-of-view _(rad)_ of the Earth's surface
  /// an observer would have at this location.
  double fieldOfView() =>
      acos(Earth.radiusMean / (Earth.radiusMean + altitude));

  /// Return `true` if this location is visible to an observer at the provided
  /// [Geodetic] coordinates.
  bool sight(final Geodetic g,
      {final AngularDistanceMethod method = AngularDistanceMethod.haversine}) {
    final fov = max(fieldOfView(), g.fieldOfView());
    return angle(g, method: method) <= fov;
  }
}
