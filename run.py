import uvicorn

if __name__ == "__main__":
    uvicorn.run(
        "app.server:APP",
        host="0.0.0.0",
        port=9096,
        reload=True,
        reload_dirs=["app"],
        log_config=None,
    )
