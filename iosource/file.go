package iosource

import "os"

// File is a [ByteSource] backed by a local file. Its ReadAt delegates to
// (*os.File).ReadAt, which is safe for concurrent use.
type File struct {
	f    *os.File
	size int64
}

// OpenFile opens a local file as a [ByteSource].
func OpenFile(path string) (*File, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	fi, err := f.Stat()
	if err != nil {
		f.Close()
		return nil, err
	}
	return &File{f: f, size: fi.Size()}, nil
}

// NewFile wraps an already-open *os.File as a [ByteSource]. The returned File
// takes ownership of f and closes it on Close.
func NewFile(f *os.File) (*File, error) {
	fi, err := f.Stat()
	if err != nil {
		return nil, err
	}
	return &File{f: f, size: fi.Size()}, nil
}

// ReadAt implements io.ReaderAt.
func (f *File) ReadAt(p []byte, off int64) (int, error) { return f.f.ReadAt(p, off) }

// Size reports the file length recorded at open time.
func (f *File) Size() (int64, error) { return f.size, nil }

// Close closes the underlying file.
func (f *File) Close() error { return f.f.Close() }
