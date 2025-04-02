import os
import sqlite3
import numpy as np
from pathlib import Path

from .constants import Constants


class SQLiteManager:
    def __init__(self, project_name: str, base_dir: Path = Constants.path.DATABASES_PATH) -> None:
        # Ensure the base directory exists
        os.makedirs(base_dir, exist_ok=True)

        # Construct the database path
        db_path: Path = base_dir / f"{project_name}.db"

        self.conn: sqlite3.Connection = sqlite3.connect(db_path)
        self._create_table()

    def _create_table(self) -> None:
        with self.conn:
            self.conn.execute("""
                CREATE TABLE IF NOT EXISTS coordinates (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    project TEXT NOT NULL,
                    element TEXT NOT NULL,
                    structure_folder TEXT,
                    structure_type TEXT,
                    idx INTEGER,
                    coordinates BLOB NOT NULL
                )
            """)

    def insert_coordinates(
        self,
        project: str,
        element: str,
        coordinates: np.ndarray,
        structure_folder: str | None = None,
        structure_type: str | None = None,
        index: int | None = None
    ) -> None:
        coords_blob: bytes = coordinates.tobytes()
        with self.conn:
            self.conn.execute("""
                INSERT INTO coordinates (project, element, structure_folder, structure_type, idx, coordinates)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (project, element, structure_folder, structure_type, index, coords_blob))

    def fetch_coordinates(
        self,
        project: str | None = None,
        element: str | None = None
    ) -> list[dict]:
        query: str = "SELECT * FROM coordinates"
        conditions: list[str] = []
        values: list[str] = []

        if project:
            conditions.append("project = ?")
            values.append(project)
        if element:
            conditions.append("element = ?")
            values.append(element)

        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        cursor: sqlite3.Cursor = self.conn.execute(query, values)
        records: list[tuple[str, str, str, str, str, int, bytes]] = cursor.fetchall()

        result: list[dict] = []
        for row in records:
            _, project, element, folder, stype, idx, coords_blob = row
            coords: np.ndarray = np.frombuffer(coords_blob, dtype=np.float64).reshape(-1, 3)
            result.append({
                "project": project,
                "element": element,
                "structure_folder": folder,
                "structure_type": stype,
                "index": idx,
                "coordinates": coords,
            })

        return result

    def close(self) -> None:
        self.conn.close()
