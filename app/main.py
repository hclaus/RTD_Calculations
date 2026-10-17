from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from app.api import routes_chemical, routes_meta, routes_microbiology

STATIC_DIR = Path(__file__).parent.parent / "static"

app = FastAPI(title="RTD Reactor Calculator")

app.include_router(routes_microbiology.router, prefix="/api/microbiology")
app.include_router(routes_chemical.router, prefix="/api/chemical")
app.include_router(routes_meta.router, prefix="/api")

app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


@app.get("/")
def index():
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/help")
def help_page():
    return FileResponse(STATIC_DIR / "help.html")
