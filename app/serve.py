#!/usr/bin/env python3
"""Serve the repository so the browser app can fetch the package source.

    python app/serve.py          then open http://localhost:8000/app/

The app has to be served over HTTP rather than opened as a file:// URL,
because it fetches swells/*.py at runtime and writes them into Pyodide's
filesystem. That is what keeps one copy of the physics instead of two.
"""

import functools
import http.server
import os
import socketserver
import webbrowser

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PORT = int(os.environ.get("PORT", 8000))


class Handler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Never cache: editing a module and hitting reload should show the change.
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def log_message(self, fmt, *args):
        pass


if __name__ == "__main__":
    socketserver.TCPServer.allow_reuse_address = True
    handler = functools.partial(Handler, directory=ROOT)
    with socketserver.TCPServer(("", PORT), handler) as httpd:
        url = f"http://localhost:{PORT}/app/"
        print(f"serving {ROOT}\n  ->  {url}\nCtrl-C to stop")
        try:
            webbrowser.open(url)
        except Exception:
            pass
        httpd.serve_forever()
