import 'package:pious_squid/src/time/epoch.dart';

/// International Atomic Time (TAI) epoch.
class EpochTAI extends Epoch {
  /// Create a new [EpochTAI] object.
  EpochTAI(final double posix) : super(posix);
}

/// Terrestrial Time (TT) epoch.
class EpochTT extends Epoch {
  /// Create a new [EpochTT] object.
  EpochTT(final double posix) : super(posix);
}

/// Barycentric Dynamical Time (TDB) epoch.
class EpochTDB extends Epoch {
  /// Create a new [EpochTDB] object.
  EpochTDB(final double posix) : super(posix);
}
