from typing import Optional

from pydantic import BaseSettings, AnyUrl


class Settings(BaseSettings):
    openapi_server_url: AnyUrl = "http://localhost:9096"
    openapi_server_maturity: str = "development"
    openapi_server_location: str = "RENCI"
    trapi_version: str = "1.5.0"
    # Redis is still used by the novelty scorer (publication lookups, db 1).
    redis_host: str = "localhost"
    redis_port: int = 6380
    redis_password: str = "supersecretpassword"
    # Clinical evidence edges are served from a memory-mapped LMDB store.
    lmdb_path: str = "./data/clinical_evidence.mdb"

    jaeger_enabled: bool = False
    jaeger_host: str = "jaeger"
    jaeger_port: int = 6831
    otel_service_name: str = "ANSWER-APPRAISER"
    otel_use_console_exporter: bool = False

    class Config:
        env_file = ".env"


settings = Settings()
