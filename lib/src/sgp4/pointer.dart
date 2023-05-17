/// Makeshift pointer object for SGP4.
class Pointer<T> {
  /// Create a [Pointer] containing the provided [value].
  Pointer(this.value);

  /// The pointer value.
  T value;
}
