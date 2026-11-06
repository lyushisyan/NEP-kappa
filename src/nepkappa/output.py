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
        elapsed = format_activity_elapsed(seconds)
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


def format_duration(seconds):
    """Format elapsed seconds as a compact human-readable duration."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = seconds % 60
    if hours:
        return f"{hours}h {minutes}m {secs:.2f}s"
    if minutes:
        return f"{minutes}m {secs:.2f}s"
    return f"{secs:.2f}s"


def format_activity_elapsed(seconds):
    """Format elapsed seconds like tqdm's compact clock."""
    total = int(seconds)
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{secs:02d}"
    return f"{minutes:02d}:{secs:02d}"
