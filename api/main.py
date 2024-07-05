from fastapi import FastAPI
from starlette.middleware.cors import CORSMiddleware

from backend.config import REACT_JS_REQUEST_ORIGIN
from api.playground_router import router as playground_router
from api.solver_router import router as solver_router
from api.classifier_router import router as classifier_router
from api.general_router import router as general_router

app = FastAPI()

# Add middleware to interface with React frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=[REACT_JS_REQUEST_ORIGIN],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
def health_check() -> str:
    return "Health check: Server is running."


app.include_router(general_router)
app.include_router(playground_router)
app.include_router(solver_router)
app.include_router(classifier_router)
