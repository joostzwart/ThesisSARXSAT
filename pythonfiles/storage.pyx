import sqlite3 as lite
import logging
import numpy as np
import io

def adapt_array(arr):
    """
    Stores numpy arrays in a database
    """
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return lite.Binary(out.read())

def convert_array(text):
    """
    Retrieves numpy arrays from a database
    """
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)

cpdef datastorage(models,parameters,name,description,datafiterror,mse,switchingsequence,time):
    """Stores the gathered data in a database
    
    Arguments:
        models {list} -- List of models.
        parameters {object} -- Parameters used during identification.
        name {string} -- Name of the experiment.
        description {string} -- Description of the experiment.
    """

    ## Register numpy arrays as a datatype
    lite.register_adapter(np.ndarray, adapt_array)
    lite.register_converter("array", convert_array)

    ## Try to connect to the database
    try:
        con = lite.connect('DataOld.db',detect_types=lite.PARSE_DECLTYPES)
    except:
        logging.error("cannot connect to database")
        exit()
    cur = con.cursor()
    cur.execute("PRAGMA foreign_keys = ON;")

    ## Insert information about the run into the table Runs
    cur.execute("INSERT INTO Runs (dataset,description,datafiterror,mse,time) values (?,?,?,?,?)",
                (name,description,datafiterror,mse,time))

    ## Identify the number of succesfull runs already in the database
    cur.execute("SELECT COUNT(*) FROM RUNS")
    run_id=cur.fetchone()

    ## Add identified models to the Models table
    for m in models:
        cur.execute("INSERT INTO Models (models,modelid) values (?,?)",
                    (np.array(m),run_id[0]))

    ## Add model information to the table Modelinfo
    cur.execute("INSERT INTO Modelinfo (id,switchingsequence,nu,ny,inputs,outputs,dw,delta) values (?,?,?,?,?,?,?,?)",
                (run_id[0],np.array(switchingsequence),parameters.nu,parameters.ny,parameters.inputs,
                 parameters.outputs,parameters.dw,parameters.delta))

    ## commit all changes and close the database
    con.commit()
    con.close()

def action(action):
    """
    This function can be used to delete all databases, create new ones or retrieve information from existing databases.
    """
    ## Register numpy arrays as a datatype
    lite.register_adapter(np.ndarray, adapt_array)
    lite.register_converter("array", convert_array)
    ## Try to connect to the database
    try:
        con = lite.connect('DataOld.db',detect_types=lite.PARSE_DECLTYPES)
    except:
        logging.error("cannot connect to database")
        exit()

    data=0
    cur = con.cursor()
    if action=="delete":
        ## Delete tables (only for testing)
        cur.execute("DROP TABLE IF EXISTS Models")
        cur.execute("DROP TABLE IF EXISTS Modelinfo")
        cur.execute("DROP TABLE IF EXISTS Runs")

    elif action=="create":
        cur.execute("CREATE TABLE Runs(id INTEGER PRIMARY KEY, dataset TEXT, description Text,"
                " datafiterror REAL,mse REAL,time REAL, date DEFAULT CURRENT_TIMESTAMP )")
        cur.execute("CREATE TABLE Models (modelid INTEGER, models array,FOREIGN KEY (modelid) REFERENCES Runs(id))")
        cur.execute("CREATE TABLE modelinfo (id INTEGER,switchingsequence array,nu INTEGER,ny INTEGER,inputs INTEGER,"
                "outputs INTEGER,dw INTEGER,delta REAL,FOREIGN KEY (id) REFERENCES Runs(id))")

    elif action=="select":
        cur.execute("SELECT * FROM Runs")
        data=cur.fetchall()
        print(data)

    elif action=="selectmodel":
        cur.execute("SELECT * FROM Models")
        data=cur.fetchall()
        print(data)

    else:
        print("Action is not recognized. Try again")
    con.commit()
    con.close()
    return data