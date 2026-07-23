from typing import Optional

from pydantic import BaseSettings, AnyUrl


class Settings(BaseSettings):
    openapi_server_url: AnyUrl = "http://localhost:9096"
    openapi_server_maturity: str = "development"
    openapi_server_location: str = "RENCI"
    trapi_version: str = "1.5.0"
    # Static lookup tables are served from memory-mapped LMDB stores.
    # Clinical evidence edges keyed by "{subject}_{object}".
    lmdb_path: str = "./data/clinical_evidence.mdb"
    # Publication years keyed by publication id (used by the novelty scorer).
    publications_lmdb_path: str = "./data/publications.mdb"

    jaeger_enabled: bool = False
    jaeger_host: str = "jaeger"
    jaeger_port: int = 6831
    otel_service_name: str = "ANSWER-APPRAISER"
    otel_use_console_exporter: bool = False

    class Config:
        env_file = ".env"


settings = Settings()
