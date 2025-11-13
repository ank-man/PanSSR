# panssrator/database.py
import sqlite3
from panssrator import config, utils

def init_db(db_path: str = config.DATABASE_FILE) -> sqlite3.Connection:
    """Initialize (or connect to) the SQLite database."""
    conn = sqlite3.connect(db_path)
    return conn

def create_marker_table(conn: sqlite3.Connection):
    """Create a table for storing markers if it does not exist."""
    sql = """
    CREATE TABLE IF NOT EXISTS markers (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        chrom TEXT,
        start INTEGER,
        end INTEGER,
        motif TEXT,
        repeat_count INTEGER,
        annotation TEXT,
        primer_info TEXT,
        epcr_results TEXT
    );
    """
    conn.execute(sql)
    conn.commit()

def insert_marker(conn: sqlite3.Connection, marker: dict):
    """Insert a marker record into the database."""
    sql = """
    INSERT INTO markers (chrom, start, end, motif, repeat_count, annotation, primer_info, epcr_results)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?);
    """
    data = (
        marker.get("chrom", ""),
        marker.get("start", 0),
        marker.get("end", 0),
        marker.get("motif", ""),
        marker.get("repeat_count", 0),
        str(marker.get("annotation", "")),
        str(marker.get("primers", "")),
        str(marker.get("amplicon_sizes", ""))
    )
    conn.execute(sql, data)
    conn.commit()

if __name__ == '__main__':
    conn = init_db()
    create_marker_table(conn)
    # Insert a test marker record
    test_marker = {
        "chrom": "chr1",
        "start": 100,
        "end": 140,
        "motif": "AT",
        "repeat_count": 10,
        "annotation": "intergenic",
        "primers": {"left": "ATGC", "right": "GCAT"},
        "amplicon_sizes": [200]
    }
    insert_marker(conn, test_marker)
    utils.logger.info("Marker inserted successfully.")
    conn.close()

