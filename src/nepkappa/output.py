"""Output details retained in the CLI log, or shown with --verbose."""

from contextlib import contextmanager
import sys
import threading
import time


def print_detail(*args, **kwargs):
    """Print normally outside the CLI; use its detail stream when available."""
    kwargs.setdefault("file", getattr(sys.stdout, "detail_stream", sys.stdout))
    print(*args, **kwargs)


@contextmanager
def activity_status(label, *, enabled=True):
    """Refresh one terminal-only elapsed-time line while the caller works."""
    stream = getattr(sys.stderr, "terminal_stream", sys.stderr)
    if not enabled or not stream.isatty():
        yield
        return

    started = time.monotonic()
    stopped = threading.Event()
    width = 0

    def draw():
        nonlocal width
        seconds = int(time.monotonic() - started)
        minutes, seconds = divmod(seconds, 60)
        hours, minutes = divmod(minutes, 60)
        elapsed = f"{hours}:{minutes:02d}:{seconds:02d}" if hours else f"{minutes:02d}:{seconds:02d}"
        text = f"{label}: running ({elapsed} elapsed)"
        width = max(width, len(text))
        stream.write("\r" + text.ljust(width))
        stream.flush()

    def refresh():
        while not stopped.wait(1.0):
            draw()

    draw()
    ticker = threading.Thread(target=refresh, daemon=True)
    ticker.start()
    try:
        yield
    finally:
        stopped.set()
        ticker.join()
        stream.write("\r" + " " * width + "\r")
        stream.flush()
