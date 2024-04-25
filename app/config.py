from typing import Optional

from pydantic import BaseSettings, AnyUrl


class Settings(BaseSettings):
    openapi_server_url: AnyUrl = "http://localhost:9096"
    openapi_server_maturity: str = "development"
    openapi_server_location: str = "RENCI"
    trapi_version: str = "1.5.0"
    redis_host: str = "localhost"
    redis_port: int = 6379
    redis_password: str = "supersecretpassword"

    class Config:
        env_file = ".env"


settings = Settings()
