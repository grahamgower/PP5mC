
import os, os.path
import gzip
from subprocess import Popen, PIPE

class GZError(Exception):
    pass

class gzopen:
    """
    Pipe file through a gzip subprocess.

    This is substantially faster than just using the gzip library,
    particularly when dealing with multiple files at a time.
    """
    def __init__(self, fn, mode="r", gz=False):

        # check for existence/permission
        if ("r" in mode and not os.path.exists(fn)) or \
                ("w" in mode and not os.path.exists(os.path.dirname(os.path.dirname(fn)))):
            raise OSError(os.errno.ENOENT, os.strerror(os.errno.ENOENT), fn)

        if ("r" in mode and not os.access(fn, os.R_OK)) or \
                ("w" in mode and not os.access(os.path.dirname(os.path.dirname(fn)), os.W_OK)):
            raise OSError(os.errno.EACCES, os.strerror(os.errno.EACCES), fn)

        if not gz and not fn.endswith(".gz"):
            # uncompressed file
            self.f = open(fn, mode)
            self.close = self.f.close
        else:
            # gz file requested
            if "r" in mode:
                try:
                    self.pipe = Popen(["gzip", "-dc", fn], bufsize=-1,
                            stdout=PIPE)
                    self.f = self.pipe.stdout
                except OSError as e:
                    if e.errno == os.errno.ENOENT:
                        # no gzip command; do gzip in current process (slow)
                        self.f = gzip.open(fn, mode)
                        self.close = self.f.close
                    else:
                        raise
            elif "w" in mode:
                f_out = open(fn, mode)
                try:
                    self.pipe = Popen(["gzip", "-c"], bufsize=-1,
                            stdin=PIPE, stdout=f_out)
                    self.f = self.pipe.stdin
                except OSError as e:
                    f_out.close()
                    if e.errno == os.errno.ENOENT:
                        # no gzip command; do gzip in current process (slow)
                        self.f = gzip.open(fn, mode)
                        self.close = self.f.close
                    else:
                        raise

        self.read = self.f.read
        self.write = self.f.write

    def __enter__(self):
        return self.f

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.pipe.stdin is not None:
            self.pipe.stdin.close()
        if self.pipe.stdout is not None:
            self.pipe.stdout.close()
        self.pipe.wait()
        if self.pipe.returncode != 0:
            raise GZError("gzip exited with return code {}".format(self.pipe.returncode))
