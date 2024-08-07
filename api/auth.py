from fastapi.security import APIKeyHeader
from fastapi import HTTPException, Security
from fastapi import Depends
import os

api_key_header = APIKeyHeader(name="x-api-Key")


def get_correct_api_key() -> str:
    api_key = os.getenv("REACT_APP_API_KEY")
    if api_key is None:
        raise HTTPException(
            status_code=500, detail="API_KEY not set in environment variables"
        )
    return api_key


def verify_api_key(api_key: str = Security(api_key_header)) -> None:
    correct_api_key = get_correct_api_key()
    if api_key != correct_api_key:
        raise HTTPException(status_code=401, detail="Invalid or missing API key.")


auth_dependencies = (
    [Depends(verify_api_key)] if os.getenv("REACT_APP_IS_CLOUD_DEPLOYMENT") == "true" else None
)
