from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional, Tuple

FrameStatus = Literal["ok", "warning95", "warning99"]
OutputFormat = Literal["ps", "pdf"]


@dataclass(frozen=True)
class ChainSummary:
    chain_id: str
    start_residue: int
    end_residue: int


@dataclass(frozen=True)
class FrameScore:
    chain_id: str
    center_residue: int
    error_value: float
    status: FrameStatus


@dataclass(frozen=True)
class ReportPaths:
    log_path: Path
    plot_path: Path
    output_format: OutputFormat


@dataclass(frozen=True)
class AnalysisResult:
    protein_id: str
    input_path: Path
    scored_frame_count: int
    rejected_frame_count: int
    rejected_frame_ratio: Optional[float]
    overall_quality_factor: Optional[float]
    average_probability: Optional[float]
    below_interaction_limit_frames: Tuple[int, ...]
    chain_summaries: Tuple[ChainSummary, ...]
    frame_scores: Tuple[FrameScore, ...]
    messages: Tuple[str, ...]
    log_text: str
    report_paths: Optional[ReportPaths] = None

    @property
    def flagged_frames(self) -> Tuple[FrameScore, ...]:
        return tuple(frame for frame in self.frame_scores if frame.status != "ok")
