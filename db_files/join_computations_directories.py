# this script takes all computations from database files and sends them to one flat directory 
import torinanet as tn
import os
import sqlite3

BASE_HASH_FUNCTION = tn.core.HashedCollection.RdkitFingerprintGenerator(tn.core.HashedCollection.FingerprintAlgorithms.RDKIT, 256, 4)

def specie_hash_func(specie):
    return str(hex(BASE_HASH_FUNCTION(specie)))

def smiles_in_db(db_session, smiles) -> bool:
    """Find the ID of a specie in a databsae session"""
    print(smiles)
    specie = tn.core.Specie(smiles)
    try:
        hash_key = specie_hash_func(specie)
    except:
        print("errors parsing {}. probably its has illegal valence".format(smiles))
        # returns true to avoid adding it to the DB
        return True
    # checks by hash_key if specie exists in the specie table
    res_smiles_strings = [s[0] for s in db_session.execute("SELECT smiles FROM species WHERE hash_key=\"{}\"".format(hash_key)).fetchall()]
    if len(res_smiles_strings) == 0:
        return False
    # first, we ty and see if there is an exact string match
    for res_smiles in res_smiles_strings:
        if smiles == res_smiles:
            return True
    # if no exact string match was found, we go to full graph comparision    
    specie_ac = tn.core.BinaryAcMatrix.from_specie(specie)
    for res_smiles in res_smiles_strings:
        ac = tn.core.BinaryAcMatrix.from_specie(tn.core.Specie(res_smiles))
        if ac == specie_ac:
            return True
    return False

def format_sql_value(x):
    if type(x) is str:
        return "\"{}\"".format(x)
    elif x is None:
        return "null"
    else:
        return str(x)

def smiles_insert_command(db_session, sid: int) -> str:
    """return the SQL insert command to insert a single specie"""
    values = [format_sql_value(s) for s in db_session.execute("SELECT * FROM species WHERE \"id\"={}".format(sid)).fetchall()[0]][1:]
    columns = [s[0] for s in db_session.execute("SELECT name FROM PRAGMA_TABLE_INFO(\"species\")").fetchall()][1:]
    return "INSERT INTO species ({}) VALUES ({})".format(",".join(columns), ",".join(values)) 

def move_computaion_files(xyz_file: str, in_file: str, out_file: str, comp_dir: str):
    """on a given db, change the location of the specie's XYZ, IN and OUT files. returns new paths to XYZ, IN and OUT files"""
    # extracting relevant information
    name = os.path.split(xyz_file)[-1][:-4]
    out_dir = os.path.split(out_file)[0]
    target_dir = os.path.join(comp_dir, name)
    if not os.path.isdir(target_dir):
        os.mkdir(target_dir)
    # running some shell commands to copy everything
    os.system("cp -r {}/* {}".format(out_dir, target_dir))
    os.system("cp {} {}".format(xyz_file, target_dir))
    os.system("cp {} {}".format(in_file, target_dir))
    return os.path.join(target_dir, name + '.xyz'), os.path.join(target_dir, name + '.inp'), os.path.join(target_dir, name + '.out')

def update_file_entries(db_session, sid: int, xyz_file: str, in_file: str, out_file: str):
    db_session.execute("UPDATE species SET xyz={}, comp_input={}, comp_output={} WHERE \"id\"={}".format(xyz_file, in_file, out_file, sid))


if __name__ == "__main__":
    db_files = ["nh3+o2.db", "nh3+o2+h2.db", "nh3+o2+ch4.db", "ch3ch3.db", "ch4+h2o.db"]
    joined_db = "joined.db"
    comp_dir = "../computation_files"
    if os.path.isfile(joined_db):
        os.remove(joined_db)
    # fetching specie table schema from the first DB file
    con = sqlite3.connect(db_files[0])
    schema = con.execute("SELECT sql FROM sqlite_schema WHERE name=\"species\"").fetchone()[0]
    con.close()
    # creating new database + specie table
    joined_con = sqlite3.connect(joined_db)
    joined_con.execute(schema)
    # now joining all DB files 
    for db_file in db_files:
        print("JOINING WITH", db_file)
        con = sqlite3.connect(db_file)
        for sid, smiles in con.execute("SELECT \"id\", smiles FROM species").fetchall():
            if not smiles_in_db(joined_con, smiles):
                query = smiles_insert_command(con, sid)
                print(query)
                joined_con.execute(query)
        con.close()
    joined_con.commit()
    # and now, moving all files to joined directory
    for sid, xyz, infile, outfile in joined_con.execute("SELECT \"id\", xyz, comp_input, comp_output FROM species"):
        print("copying specie number", sid)
        new_xyz, new_infie, new_outfie = move_computaion_files(xyz, infile, outfile, comp_dir)
        print("UPDATE species SET xyz=\"{}\", comp_input=\"{}\", comp_output=\"{}\" WHERE \"id\"={}".format(new_xyz, new_infie, new_outfie, int(sid)))
        joined_con.execute("UPDATE species SET xyz=\"{}\", comp_input=\"{}\", comp_output=\"{}\" WHERE \"id\"={}".format(new_xyz, new_infie, new_outfie, int(sid)))
    joined_con.commit()
    print("ALL DONE !")