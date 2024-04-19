/// Pseudo-pointer object.
class Pointer<T> {
  /// Create a [Pointer] containing the provided [value].
  Pointer(this.value);

  /// The pointer value.
  T value;

  @override
  String toString() => value.toString();
}
