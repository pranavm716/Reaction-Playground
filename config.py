import json
from dataclasses import dataclass


@dataclass
class Config:
    all_reactions_file_path: str
    multi_step_react_mode: bool
    max_num_solver_steps: int


def read_config(config_file: str) -> Config:
    with open(config_file) as f:
        data = json.load(f)
        return Config(**data)
