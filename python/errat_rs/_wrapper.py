from __future__ import annotations

from os import PathLike, fspath
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union, cast

from . import _native
from .models import (
    AnalysisResult,
    ChainSummary,
    FrameScore,
    FrameStatus,
    OutputFormat,
    ReportPaths,
)

Pathish = Union[str, PathLike[str]]


def analyze(
    input_path: Pathish,
    *,
    protein_id: Optional[str] = None,
    use_mmap: bool = False,
) -> AnalysisResult:
    payload = cast(
        Dict[str, Any],
        _native.analyze(
            fspath(input_path),
            protein_id=protein_id,
            use_mmap=use_mmap,
        ),
    )
    return _analysis_from_payload(payload)


def write_report(
    input_path: Pathish,
    output_dir: Pathish,
    *,
    protein_id: Optional[str] = None,
    output_format: OutputFormat = "ps",
    use_mmap: bool = False,
) -> ReportPaths:
    normalized = _normalize_output_format(output_format)
    payload = cast(
        Dict[str, Any],
        _native.write_report(
            fspath(input_path),
            fspath(output_dir),
            protein_id=protein_id,
            output_format=normalized,
            use_mmap=use_mmap,
        ),
    )
    return ReportPaths(
        log_path=Path(payload["log_path"]),
        plot_path=Path(payload["plot_path"]),
        output_format=cast(OutputFormat, payload["output_format"]),
    )


def analyze_and_write(
    input_path: Pathish,
    output_dir: Pathish,
    *,
    protein_id: Optional[str] = None,
    output_format: OutputFormat = "ps",
    use_mmap: bool = False,
) -> AnalysisResult:
    normalized = _normalize_output_format(output_format)
    payload = cast(
        Dict[str, Any],
        _native.analyze_and_write(
            fspath(input_path),
            fspath(output_dir),
            protein_id=protein_id,
            output_format=normalized,
            use_mmap=use_mmap,
        ),
    )
    return _analysis_from_payload(payload)


def _analysis_from_payload(payload: Dict[str, Any]) -> AnalysisResult:
    chain_summaries = tuple(
        ChainSummary(
            chain_id=str(item["chain_id"]),
            start_residue=int(item["start_residue"]),
            end_residue=int(item["end_residue"]),
        )
        for item in cast(Iterable[Dict[str, Any]], payload["chain_summaries"])
    )
    frame_scores = tuple(
        FrameScore(
            chain_id=str(item["chain_id"]),
            center_residue=int(item["center_residue"]),
            error_value=float(item["error_value"]),
            status=cast(FrameStatus, item["status"]),
        )
        for item in cast(Iterable[Dict[str, Any]], payload["frame_scores"])
    )
    return AnalysisResult(
        protein_id=str(payload["protein_id"]),
        input_path=Path(payload["input_path"]),
        scored_frame_count=int(payload["scored_frame_count"]),
        rejected_frame_count=int(payload["rejected_frame_count"]),
        rejected_frame_ratio=_optional_float(payload["rejected_frame_ratio"]),
        overall_quality_factor=_optional_float(payload["overall_quality_factor"]),
        average_probability=_optional_float(payload["average_probability"]),
        below_interaction_limit_frames=tuple(
            int(frame)
            for frame in cast(Iterable[Any], payload["below_interaction_limit_frames"])
        ),
        chain_summaries=chain_summaries,
        frame_scores=frame_scores,
        messages=tuple(
            str(message) for message in cast(Iterable[Any], payload["messages"])
        ),
        log_text=str(payload["log_text"]),
        report_paths=_optional_report_paths(payload.get("report_paths")),
    )


def _optional_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    return float(value)


def _normalize_output_format(output_format: str) -> OutputFormat:
    normalized = output_format.lower()
    if normalized not in {"ps", "pdf"}:
        raise ValueError("output_format must be either 'ps' or 'pdf'")
    return cast(OutputFormat, normalized)


def _optional_report_paths(value: Any) -> Optional[ReportPaths]:
    if value is None:
        return None
    payload = cast(Dict[str, Any], value)
    return ReportPaths(
        log_path=Path(payload["log_path"]),
        plot_path=Path(payload["plot_path"]),
        output_format=cast(OutputFormat, payload["output_format"]),
    )
