from fastapi import APIRouter
from backend.reaction import Reaction
from backend.computations import ALL_REACTIONS
from api.auth import auth_dependencies


router = APIRouter(tags=["general"], dependencies=auth_dependencies)


@router.get("/reactions")
def all_reactions() -> list[Reaction]:
    """Returns a list of all available reactions."""
    return list(ALL_REACTIONS.values())
