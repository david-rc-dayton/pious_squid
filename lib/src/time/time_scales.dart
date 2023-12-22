import 'package:pious_squid/src/time/epoch.dart';

/// International Atomic Time (TAI) epoch.
class EpochTAI extends Epoch {
  /// Create a new [EpochTAI] object.
  EpochTAI(super.posix);
}

/// Terrestrial Time (TT) epoch.
class EpochTT extends Epoch {
  /// Create a new [EpochTT] object.
  EpochTT(super.posix);
}

/// UT1 epoch.
class EpochUT1 extends Epoch {
  /// Create a new [EpochUT1] object.
  EpochUT1(super.posix);
}

/// Barycentric Dynamical Time (TDB) epoch.
class EpochTDB extends Epoch {
  /// Create a new [EpochTDB] object.
  EpochTDB(super.posix);
}
