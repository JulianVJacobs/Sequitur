import sqlite3

def merge_tables():
    conn = sqlite3.connect('/workspace/sequitur/data/input/ERS218597/suffix.db')
    cursor = conn.cursor()
    
    # Count the number of records in both tables
    cursor.execute('SELECT COUNT(*) FROM ERS218597')
    count_ERS218597 = cursor.fetchone()[0]
    
    cursor.execute('SELECT COUNT(*) FROM ERR234359')
    count_ERR234359 = cursor.fetchone()[0]
    
    # Determine which table to drop and which to rename
    if count_ERS218597 >= count_ERR234359:
        # Drop ERR234359 and keep ERS218597
        cursor.execute('DROP TABLE IF EXISTS ERR234359')
    else:
        # Drop ERS218597 and rename ERR234359 to ERS218597
        cursor.execute('DROP TABLE IF EXISTS ERS218597')
        cursor.execute('ALTER TABLE ERR234359 RENAME TO ERS218597')
    
    # Commit the changes
    conn.commit()
    
    # Vacuum the database to free up space
    cursor.execute('VACUUM')
    
    # Close the connection
    conn.close()

if __name__ == "__main__":
    merge_tables()