"""OpenTelemetry tracing setup for the FastAPI app.

Must be called *after* the gunicorn fork (from a post_fork hook) so the
BatchSpanProcessor background thread and gRPC channel are created in the
worker process, not the master. See
https://oneuptime.com/blog/post/2026-02-06-troubleshoot-fastapi-uvicorn-reload/view
"""


def setup_otel(worker_pid=None):
    from .config import settings

    if settings.jaeger_enabled in ("false", "False"):
        return

    from opentelemetry import trace
    from opentelemetry.sdk.resources import SERVICE_NAME, Resource
    from opentelemetry.sdk.trace import TracerProvider
    from opentelemetry.sdk.trace.export import BatchSpanProcessor
    from opentelemetry.instrumentation.fastapi import FastAPIInstrumentor
    from opentelemetry.instrumentation.httpx import HTTPXClientInstrumentor

    resource_attrs = {SERVICE_NAME: settings.otel_service_name}
    if worker_pid is not None:
        resource_attrs["worker.pid"] = worker_pid
    resource = Resource.create(attributes=resource_attrs)
    provider = TracerProvider(resource=resource)

    if settings.otel_use_console_exporter:
        from opentelemetry.sdk.trace.export import ConsoleSpanExporter
        processor = BatchSpanProcessor(ConsoleSpanExporter())
    else:
        from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import OTLPSpanExporter
        otlp_host = settings.jaeger_host.rstrip('/')
        otlp_port = settings.jaeger_port
        processor = BatchSpanProcessor(OTLPSpanExporter(endpoint=f'{otlp_host}:{otlp_port}'))

    provider.add_span_processor(processor)
    trace.set_tracer_provider(provider)

    from .server import APP
    FastAPIInstrumentor.instrument_app(APP, excluded_urls="docs,openapi.json,redis_ready")
    HTTPXClientInstrumentor().instrument()
