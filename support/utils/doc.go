// Package utils provides small, general-purpose helpers used across the
// toolkit.
//
// It includes:
//
//   - [Semaphore], a counting semaphore backed by a buffered channel for
//     bounding concurrency.
//   - [PositionTrackingReader], an [io.Reader] wrapper that counts the bytes
//     read through it.
//   - [TrimFloat], which formats a float and strips trailing zeros.
package utils
