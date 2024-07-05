from fastapi import APIRouter, Depends
from backend.reaction import Reaction
from backend.computations import ALL_REACTIONS
from api.auth import verify_api_key


router = APIRouter(tags=["general"], dependencies=[Depends(verify_api_key)])


@router.get("/reactions")
def all_reactions() -> list[Reaction]:
    """Returns a list of all available reactions."""
    return list(ALL_REACTIONS.values())
