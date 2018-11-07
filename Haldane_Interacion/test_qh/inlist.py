Dyson={
"Control": {
    "__Execute" : ["python", "./dyson/main.py"],
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : True,
    "__PBSCommand": "#PBS -l mem=5gb"
    },
"Job": {
    "DysonOnly": True,
    "SumRule": False 
    }
}

Beta=2
Order=3
Common={
"Tau": {"MaxTauBin" :512, "Beta": Beta},
"Lattice":  {
    #2D lattice
    #"Name": "Square_F", "NSublat": 4,
    #"Name": "Checkerboard", "NSublat": 2,
    #"Name": "ValenceBond", "NSublat": 2,
    "Name": "Haldane", "NSublat": 2,
    #"Name": "Honeycomb", "NSublat": 2,
    #"Name": "Kagome", "NSublat": 3,
    "L": [16,16]

    #3D lattice
    #"Name": "Cubic", "NSublat": 1,
    #"Name": "3DCheckerboard", "NSublat": 2,
    #"Name": "Pyrochlore", "NSublat": 4,
    #"L": [4,4,4]
    },
"Model": {
    "Name": "Hubbard",
    #"Description": ["ImW",],
    "Interaction": [1.0,0.0,0.0],
    "ExternalField": [0.00, 0.00, 0.0, 0.0],
    "ChemicalPotential": -1.00,
    "Hopping": [1.0,0.0,0.0]
    #ExternalField on Sublattice A and B
    },
}

Dyson["Dyson"]={
    "SleepTime": 90,
    #"SleepTime": 300,
    "Annealing": {
        "DeltaField": [0.0, -0.0, 0.0, 0.0],
        "Interval": [-0.1, 0.1, -0.0, -0.0]
        }
    }

import job_class as job
'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
TO_DO = []
Dyson.update(Common)
TO_DO.append(job.JobDyson(Dyson))
CPU = 4
SLEEP = 1    #check job status for every SLEEP seconds
