#include <ilcp/cp.h>
#include <ilcp/cpext.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

typedef IloArray <IloIntVarArray> IloIntVarArray2;

class FileError: public IloException {
public:
  FileError() : IloException("Cannot open data file") {}
};


// Define new global constraint
class BATCH : public IlcConstraintI {

  private:

	IlcIntervalVarArray _activityIntervals, _batchIntervals;
	IlcIntVarArray _batchDecisions, _loads;

  private:

	   // filtering algorithms
     void synchronize(IlcIntervalVar activity, IlcIntervalVar batch) const;
	   void overlap(IlcIntervalVar activity, IlcIntervalVar batch) const;
	   void tighten(IlcIntervalVar activity, IlcIntervalVar batch) const;

  public:

  	BATCH(IloCP cp, IlcIntervalVarArray activityIntervals, IlcIntervalVarArray batchIntervals, IlcIntVarArray batchDecisions, IlcIntVarArray loads);
  	~BATCH(){}; // empty destructor

    virtual void propagate ();
    virtual void post();

	  // propagation demons
    void propagateBatchAssigned(IlcIntervalVar activity, IlcIntervalVarArray batchIntervals, IlcIntVar decision);
  	void propagateTimeWindowActivityChanged(IlcIntervalVar activity, IlcIntervalVarArray batchIntervals, IlcIntVar decision);
  	void propagateTimeWindowBatchChanged(IlcInt batchIndex, IlcIntervalVar batch, IlcIntervalVarArray activityIntervals, IlcIntVarArray batchDecisions);
  	void propagateBatchesChanged(IlcIntervalVar activity, IlcIntervalVarArray batchIntervals, IlcIntVar decision);
  	void propagateLoadsFixed(IlcInt index, IlcIntVarArray loads);

};

BATCH::BATCH(IloCP cp, IlcIntervalVarArray activityIntervals, IlcIntervalVarArray batchIntervals, IlcIntVarArray batchDecisions, IlcIntVarArray loads)
	: IlcConstraintI(cp), _activityIntervals(activityIntervals), _batchIntervals(batchIntervals), _batchDecisions(batchDecisions), _loads(loads) {}

ILCCTDEMON3(Demon1,BATCH,propagateBatchAssigned,IlcIntervalVar, activity, IlcIntervalVarArray, batchIntervals, IlcIntVar, decision)
ILCCTDEMON3(Demon2,BATCH,propagateTimeWindowActivityChanged,IlcIntervalVar, activity, IlcIntervalVarArray, batchIntervals, IlcIntVar, decision)
ILCCTDEMON3(Demon3,BATCH,propagateBatchesChanged,IlcIntervalVar, activity, IlcIntervalVarArray, batchIntervals, IlcIntVar, decision)
ILCCTDEMON4(Demon4,BATCH,propagateTimeWindowBatchChanged,IlcInt, batchIndex, IlcIntervalVar, batch, IlcIntervalVarArray, activityIntervals, IlcIntVarArray, batchDecisions)
ILCCTDEMON2(Demon5,BATCH,propagateLoadsFixed,IlcInt, index, IlcIntVarArray, loads)


void BATCH::propagateBatchAssigned(IlcIntervalVar activity, IlcIntervalVarArray batchIntervals, IlcIntVar decision) {

	IlcInt assignment = decision.getValue();

	IlcStartAtStart(activity,batchIntervals[assignment],0);
	IlcEndAtEnd(activity,batchIntervals[assignment],0);

}

void BATCH::propagateTimeWindowActivityChanged(IlcIntervalVar activity, IlcIntervalVarArray batchIntervals, IlcIntVar decision) {

	// Look at only min batch and max batch in the domain of decision
	// Then if there is change in the domain, demon3 will be active
	// Then if demon3 changes the domain, activity comes back here and new min and max batch are analyzed

	if (decision.isFixed() == false) {

		IlcInt activityEST = activity.getStartMin();
		IlcInt activityLFT = activity.getEndMax();
		IlcInt activitySize = activity.getSizeMin();

		for (IlcInt i = decision.getMin(); i < decision.getMax() + 1; i++) {

			IlcIntervalVar batchInterval = _batchIntervals[i];

			IlcInt batchEST = batchInterval.getStartMin();
			IlcInt batchLFT = batchInterval.getEndMax();
			IlcInt batchSize = batchInterval.getSizeMin();

			if ((activityEST >= batchLFT) || (activityLFT <= batchEST)) {

				decision.removeValue(i);

			}

			else {


				IlcInt maxStart;
				IlcInt minEnd;
				IlcInt overlap;

				if (activityEST >= batchEST) {maxStart = activityEST;}
				else {maxStart = batchEST;}

				if (activityLFT <= batchLFT) {minEnd = activityLFT;}
				else {minEnd = batchLFT;}

				overlap = minEnd - maxStart;

				if ((overlap < activitySize) || (overlap < batchSize)) {
					decision.removeValue(i);
				}

			}

		}

	}

}

void BATCH::propagateBatchesChanged(IlcIntervalVar activity, IlcIntervalVarArray batchIntervals, IlcIntVar decision) {

  IlcIntervalVar earliestBatch = batchIntervals[decision.getMin()];
	IlcIntervalVar latestBatch = batchIntervals[decision.getMax()];

	activity.setStartMin(earliestBatch.getStartMin());
	activity.setEndMax(latestBatch.getEndMax());

}

void BATCH::propagateTimeWindowBatchChanged(IlcInt batchIndex, IlcIntervalVar batch, IlcIntervalVarArray activityIntervals, IlcIntVarArray batchDecisions) {

	IlcInt batchEST = batch.getStartMin();
	IlcInt batchLFT = batch.getEndMax();
	IlcInt batchSize = batch.getSizeMin();

	// for all activities, check if this batch is in the domain of each activity and if yes, test time windows
	for (IlcInt i = 0; i < activityIntervals.getSize(); i++) {

		IlcIntVar decision = batchDecisions[i];
		if (decision.isInDomain(batchIndex)) {

			IlcIntervalVar activity = activityIntervals[i];

			IlcInt activityEST = activity.getStartMin();
			IlcInt activityLFT = activity.getEndMax();
			IlcInt activitySize = activity.getSizeMin();

			if ((activityEST >= batchLFT) || (activityLFT <= batchEST)) {
				decision.removeValue(i);
			}

			else {
				IlcInt maxStart;
				IlcInt minEnd;
				IlcInt overlap;

				if (activityEST >= batchEST) {maxStart = activityEST;}
				else {maxStart = batchEST;}

				if (activityLFT <= batchLFT) {minEnd = activityLFT;}
				else {minEnd = batchLFT;}

				overlap = minEnd - maxStart;

				if ((overlap < activitySize) || (overlap < batchSize)) {
					decision.removeValue(i);
				}

			}

		}

	}

}

void BATCH::propagateLoadsFixed(IlcInt index, IlcIntVarArray loads) {

	// Empirically makes search slower

	if (loads[index].getValue() == 0) {
		IlcInt nbBatches = loads.getSize();

		for (IlcInt i = index + 1; i < nbBatches; i++) {
			loads[i].setMax(0);
		}

	}

}

void BATCH::propagate () {

	IlcInt nbJobs = _batchDecisions.getSize();

	for (IlcInt i = 0; i < nbJobs; i++) {
  	 IlcIntVar decision = _batchDecisions[i];

     if (decision.isFixed()) {
        IlcIntervalVar activity = _activityIntervals[i];
	  		IlcIntervalVar batch = _batchIntervals[decision.getValue()];

	  		// Try to cp.add() constraints instead of pruning it
	  		// acitivity.getCP()
	  		// Then find a way to refer to the start time of an interval
	  		// make sure this is only done once (using demons?)
	  		activity.setStartMin(batch.getStartMin());
	  		batch.setStartMin(activity.getStartMin());
	  		activity.setEndMax(batch.getEndMax());
	  		batch.setEndMax(activity.getEndMax());
	  		activity.setSizeMin(batch.getSizeMin());
	  		batch.setSizeMin(activity.getSizeMin());

	   }

	   else {
	  		IlcIntervalVar activity = _activityIntervals[i];

	  		IlcInt activityEST = activity.getStartMin();
	  		IlcInt activityLFT = activity.getEndMax();
	  		IlcInt activitySize = activity.getSizeMin();

	  		for (IlcInt j = 0; j < nbJobs; j++) {
	  			if (decision.isInDomain(j)) {
	  				IlcIntervalVar currentBatch = _batchIntervals[j];

	  				IlcInt batchEST = currentBatch.getStartMin();
	  				IlcInt batchLFT = currentBatch.getEndMax();
	  				IlcInt batchSize = currentBatch.getSizeMin();

	  				// Case 1
	  				if ((activityEST >= batchLFT) || (activityLFT <= batchEST)) {
	  					decision.removeValue(j);

	  				}

	  				// Case 2
	  				else if ((activityEST >= batchEST) && (activityLFT <= batchLFT)) {
	  					IlcInt overlap = activityLFT - activityEST;

	  					if ((overlap < activitySize) || (overlap < batchSize)) {
	  						decision.removeValue(j);

	  					}

	  				}

	  				// Case 3
	  				else if ((batchEST > activityEST) && (batchLFT < activityLFT)) {
	  					IlcInt overlap = batchLFT - batchEST;

	  					if ((overlap < activitySize) || (overlap < batchSize)) {
	  						decision.removeValue(j);

	  					}

	  				}

	  				// Case 4
	  				// Condition (batchEST < activityLFT) true since first if was false
	  				else if ((activityEST < batchEST) && (activityLFT <= batchLFT)) {
	  					IlcInt overlap = activityLFT - batchEST;

	  					if ((overlap < activitySize) || (overlap < batchSize)) {
	  						decision.removeValue(j);

	  					}

	  				}

	  				// Case 5
	  				else {
	  					IlcInt overlap = batchLFT - activityEST;

	  					if ((overlap < activitySize) || (overlap < batchSize)) {
	  						decision.removeValue(j);

	  					}

	  				}

	  			}

	  		}

	  		IlcIntervalVar earliestBatch = _batchIntervals[decision.getMin()];
	  		IlcIntervalVar latestBatch = _batchIntervals[decision.getMax()];

	  		activity.setStartMin(earliestBatch.getStartMin());
	  		activity.setEndMax(latestBatch.getEndMax());

	  		// Update LB on size of activity by setting it to smallest processing time of batches it can be assigned to

	  }

		//for all batches, if load = 0, check the load of all batches after it. if any batch has LB(load) > 0, fail

	}

}

void BATCH::post(){

	for (IlcInt i = 0; i < _batchDecisions.getSize(); i++) {

		_batchDecisions[i].whenRange(this);
		_activityIntervals[i].whenIntervalDomain(this);
		_batchIntervals[i].whenIntervalDomain(this);
		_activityIntervals[i].whenSize(this);
		_batchIntervals[i].whenSize(this);
	}

	for (IlcInt i = 0; i < _batchDecisions.getSize(); i++) {

		// When batch of activiy is assigned
		_batchDecisions[i].whenValue(Demon1(getCP(), this, _activityIntervals[i], _batchIntervals, _batchDecisions[i]));

		// When time window of activity changes
		_activityIntervals[i].whenIntervalDomain(Demon2(getCP(), this, _activityIntervals[i], _batchIntervals, _batchDecisions[i]));

		// When time window of batch changes
		_batchIntervals[i].whenIntervalDomain(Demon4(getCP(), this, i, _batchIntervals[i], _activityIntervals, _batchDecisions));

		// When domain of batches of activity changes
		_batchDecisions[i].whenRange(Demon3(getCP(), this, _activityIntervals[i], _batchIntervals, _batchDecisions[i]));

		// When load is fixed (i.e. batch is complete)
		_loads[i].whenValue(Demon5(getCP(), this, i, _loads));

	}

}

IlcConstraint Batch(IlcIntervalVarArray activityIntervals, IlcIntervalVarArray batchIntervals, IlcIntVarArray batchDecisions){

		IlcIntVar decision = batchDecisions[0];
		IloCP cp = batchDecisions[0].getCP();

    return new (cp.getHeap()) BATCH(cp,activityIntervals, batchIntervals, batchDecisions);
}


// MODEL-ENGINE CONSTRAINT WRAPPING

ILOCPCONSTRAINTWRAPPER4(IloBatch, cp,
                        IloIntervalVarArray, activities,
						            IloIntervalVarArray, batches,
                        IloIntVarArray, decisions,
						            IloIntVarArray, loads) {
  // Extraction from model data.
  // Extracting sequence, it also extract intervals of sequences
  use(cp, activities);
  use(cp, batches);
  use(cp, decisions);
  use(cp, loads);

  // Create engine data structures
  IloInt size = activities.getSize();
  IlcIntervalVarArray intervalsOfActivities(cp, size);
  IlcIntervalVarArray intervalsOfBatches(cp, size);
  IlcIntVarArray decisionsBatches(cp,size);
  IlcIntVarArray loadsBatches(cp,size);

  for(IloInt i = 0; i < size; ++i) {
    intervalsOfActivities[i] = cp.getInterval(activities[i]);
    intervalsOfBatches[i] = cp.getInterval(batches[i]);
	  decisionsBatches[i] = cp.getIntVar(decisions[i]);
	  loadsBatches[i] = cp.getIntVar(loads[i]);
  }

  // Create constraint
  return new (cp.getHeap())
    BATCH(cp, intervalsOfActivities, intervalsOfBatches, decisionsBatches, loadsBatches);
}

// End of new global constraint


int main(int argc, const char* argv[]){
  // Create environment
  IloEnv env;

  IloInt numberInstances = 48;
  IloTimer timer(env);

  timer.start();

  for (IloInt x = 0; x < numberInstances; x++) {

  for (IloInt instanceSize = 2; instanceSize < 3; instanceSize++) {

  for (IloInt trial = 0; trial < 1; trial++) {

  std::stringstream instanceFile;

  instanceFile << "../../../examples/data/jobshop/newInstance_" << x << ".txt";
  std::string data = instanceFile.str();

  std::ifstream file(data);

  std::vector<std::vector<int>> vv;
  std::string line;

  std::getline(file, line);

  std::stringstream ss(line);

  int numberJobs, numberMachines;

  ss >> numberJobs;
  ss >> numberMachines;


  while(std::getline(file, line))
  {
      std::stringstream ss(line);
      int i;
      std::vector<int> v;
      while( ss >> i )
         v.push_back(i);
      vv.push_back(v);
  }

  std::vector<std::vector<int>> processingTimesActivities;
  std::vector<std::vector<int>> machinesActivities;
  std::vector<std::vector<int>> sizesActivities;
  std::vector<int> capacitiesMachines;

  int index = 0;


  // Get processing time of jobs
  for (int i = 0; i < numberJobs; i++) {
	  std::vector<int> v1;

	  for (int j = 0; j < numberMachines; j++) {
		  v1.push_back(vv[index][j]);
	  }

	  processingTimesActivities.push_back(v1);
	  index += 1;

  }

  // Get sequence of machines for each job
  for (int i = 0; i < numberJobs; i++) {
	  std::vector<int> v1;
	  for (int j = 0; j < numberMachines; j++) {
		  v1.push_back(vv[index][j]);
	  }

	  machinesActivities.push_back(v1);
	  index += 1;

  }

  // Get sizes of jobs
  for (int i = 0; i < numberJobs; i++) {
	  std::vector<int> v1;
	  for (int j = 0; j < numberMachines; j++) {
		  v1.push_back(vv[index][j]);
	  }

	  sizesActivities.push_back(v1);
	  index += 1;

  }

  // Get capacity of each machine
  for (int i = 0; i < numberMachines; i++)
	  capacitiesMachines.push_back(vv[index][i]);

  /*
  for (int n = 0; n < instanceSize; n++) {

    for (int i = 0; i < numberJobs; i++) {

      std::vector<int> v1;
	  std::vector<int> v2;
	  std::vector<int> v3;

	  for (int j = 0; j < numberMachines; j++) {

	    v1.push_back(vv[i][2*j + 1]);

	    v2.push_back(vv[i][2*j]);

		v3.push_back(1);

	  }

	  processingTimesActivities.push_back(v1);
	  machinesActivities.push_back(v2);
	  sizesActivities.push_back(v3);
	  capacitiesMachines.push_back(instanceSize);

    }

  }*/


  printf("Number of Jobs: %d \n", numberJobs);
  printf("Number of Machines: %d \n", numberMachines);

  printf("Processing Time Job 4 Actv. 3: %d \n", processingTimesActivities[4][3]);
  printf("Machine Job 4 Actv. 3: %d \n", machinesActivities[4][3]);

  try {

	// Data file location
    const char* filename = "../../../examples/data/jobshop_default.data";

	// Fail limit
    IloInt failLimit = 1000000;

    if (argc > 1)
      filename = argv[1];
    if (argc > 2)
      failLimit = atoi(argv[2]);
    std::ifstream file(filename);
    if (!file){
      env.out() << "usage: " << argv[0] << " <file> <failLimit>" << std::endl;
      throw FileError();
    }

	// Creat model object
    IloModel model(env);

	// Number of jobs, number of machines
    IloInt nbJobs, nbMachines;
    file >> nbJobs;
    file >> nbMachines;

	nbJobs = numberJobs;
	nbMachines = numberMachines;

	IloIntVar Cmax(env);

	// Create bi-dimensional array with interval variables, of size nbMachines x nbMachines
	// We actually only need nbJobs x nbachines elements
    // IloIntervalVarArray2 machines(env, nbMachines);
	IloCumulFunctionExprArray machines(env, nbMachines);
	IloCumulFunctionExprArray packCumul(env, nbMachines);

	// For all machines
    for (IloInt j = 0; j < nbMachines; j++) {

	  // Element of machines is an interval array
      machines[j] = IloCumulFunctionExpr(env);
	  packCumul[j] = IloCumulFunctionExpr(env);

	}

	// Create interval variables for each batch. e.g. batchMatrix[1][2] is batch 2 in machine 3
	IloIntervalVarArray2 batchIntervalsMatrix(env);

	// Create integer variables for each activity = batch it is assigned in machine k
	IloIntVarArray2 batchMatrix(env);

	for (IloInt i = 0; i < nbJobs; i++) {

		IloIntervalVarArray batchIntervalsArray(env);
		IloIntVarArray batchArray(env);

		for (IloInt j = 0; j < nbMachines; j++) {

			IloIntervalVar batchInterval(env);
			batchIntervalsArray.add(batchInterval);

			IloIntVar batch(env,0,nbJobs-1);
			batchArray.add(batch);

		}

		batchIntervalsMatrix.add(batchIntervalsArray);
		batchMatrix.add(batchArray);

	}

	// Add bin packing constraint for each machine

	IloIntVarArray2 loadMatrix(env);

	for (IloInt j = 0; j < nbMachines; j++) {

		IloIntVarArray load(env);
		IloIntVarArray batchesOfMachine(env);

		// Need to get sizes of activities processed in machine j
		IloIntArray sizesOfActivities(env);

		for (IloInt i = 0; i < nbJobs; i++) {

			for (IloInt k = 0; k < nbMachines; k++) {

				IloInt machine = machinesActivities[i][k];

				if (machine == j) {

					IloInt size = sizesActivities[i][k];
					sizesOfActivities.add(size);

				}

			}

		}

		// Create variable for load of each batch in machine j
		for (IloInt i = 0; i < nbJobs; i++) {

			batchesOfMachine.add(batchMatrix[i][j]);
			IloInt maxcap = capacitiesMachines[j];
			load.add(IloIntVar(env,0,maxcap));

		}

		loadMatrix.add(load);

		// Add bin-packing global constraint
		// Load of each batch has lower bound 0 and upper bound equal to capacity of machine
		// batchesOfMachine[i] = batch of job i in machine j
		// sizesOfActivities[i] = size of job i in machine j

		model.add(IloPack(env,load,batchesOfMachine,sizesOfActivities));

		// Break simmetry: if load of a batch is zero, the load of subsequent batches must be zero
		/*for (IloInt i = 0; i < nbJobs - 1; i++) {

			for (IloInt b = i + 1; b < nbJobs; b++) {

				model.add(IloIfThen(env,load[i]==0,load[b]==0));

			}

		}*/

	}



	// Add precedence constraints between sucessive batches
	for (IloInt j = 0; j < nbMachines; j++) {

		IloIntervalVarArray batchMachine(env);

		/*for (IloInt i = 1; i < nbJobs; i++) {

			model.add(IloEndBeforeStart(env, batchIntervalsMatrix[i-1][j], batchIntervalsMatrix[i][j]));

		}*/

		for (IloInt i = 0; i < nbJobs - 1; i++) {

			batchMachine.add(batchIntervalsMatrix[i][j]);
			model.add(batchMatrix[i][j] <= i );

			for (IloInt k = i + 1; k < nbJobs; k++) {

				// incorrect symmetry breaking constraint: activity 3 could be assigned to batch 1 if activity 1 is assigned to batch 1, and activity 2 to batch 2
				//model.add(batchMatrix[i][j] <= batchMatrix[k][j]);

			}

		}

		model.add(IloNoOverlap(env,batchMachine));

		//model.add(IloEndOf(batchIntervalsMatrix[nbJobs-1][j]) == IloStartOf(batchIntervalsMatrix[nbJobs-1][j]) + IloSizeOf(batchIntervalsMatrix[nbJobs-1][j]));

	}


	// Create array of expressions to represent completion time of each job
    IloIntExprArray ends(env);

	// Create matrix to store inverval variables for each activity
	IloIntervalVarArray2 intervalsMatrix(env);

	// For all jobs
    for (IloInt i = 0; i < nbJobs; i++) {

	  IloIntervalVarArray intervalsArray(env);

	  // Create interval variable prec. If this is the first activity,
	  // it has duration 0 and no predecessors (dummy activity)
      IloIntervalVar prec;

	  // For all activities in the job
      for (IloInt j = 0; j < nbMachines; j++) {

		if (j > 0)
			intervalsArray.add(prec);

		// m is the machine of the activity, d is the duration of the next activity
        IloInt m, d;
        file >> m;
        file >> d;

		m = machinesActivities[i][j];
		d = processingTimesActivities[i][j];

		// Create interval variable for next activity in the job, with duration d
        IloIntervalVar ti(env);

		model.add(IloSizeOf(ti) >= d);

		// Adds interval variable to interval variable array machine[m], where m is the machine of the current activity
        // machines[m].add(ti);
		machines[m] += IloPulse(ti,sizesActivities[i][j]);

        if (0 != prec.getImpl())

		  // Add predecessor constraint from prec to ti
          model.add(IloEndBeforeStart(env, prec, ti));

		// Make prec equal to current activity
        prec = ti;



		// Create dummy activity that starts and ends with batch, and has size = cap - load
		IloIntervalVar dummy1(env);

		model.add(IloStartOf(dummy1) = IloStartOf(batchIntervalsMatrix[i][j]));
		model.add(IloSizeOf(dummy1) = IloSizeOf(batchIntervalsMatrix[i][j]));


		machines[j] += IloPulse(dummy1,0,capacitiesMachines[j]);
		model.add(IloHeightAtStart(dummy1,machines[j]) = capacitiesMachines[j] - loadMatrix[j][i]);

		// Creat dummy activity between batches
		IloIntervalVar dummy2(env);


		/*
		// Dummy activity before first batch
		if (i == 0) {

			//model.add(IloStartOf(dummy2) = 0);
			model.add(IloEndOf(dummy2) = IloStartOf(batchIntervalsMatrix[i][j]));
			model.add(IloSizeOf(dummy2) = IloStartOf(batchIntervalsMatrix[i][j]));

		}

		if (i > 0) {

			model.add(IloStartOf(dummy2) = IloEndOf(batchIntervalsMatrix[i-1][j]));
			model.add(IloEndOf(dummy2) = IloStartOf(batchIntervalsMatrix[i][j]));
			model.add(IloSizeOf(dummy2) = IloStartOf(batchIntervalsMatrix[i][j]) - IloEndOf(batchIntervalsMatrix[i-1][j]));

		}

		machines[j] += IloPulse(dummy2,capacitiesMachines[j]);

		// Dummy activity after last batch
		if (i == nbJobs - 1) {

			IloIntervalVar dummy2(env);

			model.add(IloStartOf(dummy2) = IloEndOf(batchIntervalsMatrix[i][j]));
			model.add(IloEndOf(dummy2) = Cmax);
			model.add(IloSizeOf(dummy2) = Cmax - IloEndOf(batchIntervalsMatrix[i][j]));

			machines[j] += IloPulse(dummy2,capacitiesMachines[j]);

		}

		// Create pack using Cumul
		IloInt processingTime = 1;

		IloIntervalVar dummy3(env,processingTime);

		model.add(IloStartOf(dummy3) = batchMatrix[i][m]);
		packCumul[m] += IloPulse(dummy3,sizesActivities[i][j]);
		*/

      }




	  intervalsArray.add(prec);

	  intervalsMatrix.add(intervalsArray);

	  // Add to expression array ends the end of the last activity
      ends.add(IloEndOf(prec));
    }

	// Add to expression array ends the end of all final batches for each machine
	for (IloInt j = 0; j < nbMachines; j++) {

		//ends.add(IloEndOf(batchIntervalsMatrix[nbJobs-1][j]));

		for (IloInt i = 0; i < nbMachines; i++) {

			ends.add(IloEndOf(batchIntervalsMatrix[i][j]));

		}

	}

	// Add logic constraint: for all machines, if an activity is assigned to a batch, they must be synchronized
	/*for (IloInt j = 0; j < nbMachines; j++) {

		for (IloInt i = 0; i < nbJobs; i++) {

			for (IloInt k = 0; k < nbJobs; k++) {


				IloInt machine = machinesActivities[i][j];

				IloAnd and(env);

				IloInt ptime = processingTimesActivities[i][j];

				and.add(IloSizeOf(batchIntervalsMatrix[k][machine]) == IloSizeOf(intervalsMatrix[i][j]));


				and.add(IloStartOf(batchIntervalsMatrix[k][machine]) == IloStartOf(intervalsMatrix[i][j]));


				model.add(IloIfThen(env,batchMatrix[i][machine] == k,and));

			}



		}

	}*/

	// For all machines
    for (IloInt j = 0; j < nbMachines; j++) {

	  // Add no overlap constraint
      model.add(machines[j] <= capacitiesMachines[j]);
	  //model.add(packCumul[j] <= capacitiesMachines[j]);

	}

	// New global constraint

	for (IloInt j = 0; j < nbMachines; j++) {

		IloIntervalVarArray activityIntervals(env);
		IloIntervalVarArray batchIntervals(env);
		IloIntVarArray batchDecisions(env);

		for (IloInt i = 0; i < nbJobs; i++) {

			batchIntervals.add(batchIntervalsMatrix[i][j]);
			batchDecisions.add(batchMatrix[i][j]);

			for (IloInt k = 0; k < nbMachines; k++) {

				if (machinesActivities[i][k] == j) {

					activityIntervals.add(intervalsMatrix[i][k]);

				}

			}

		}

		model.add(IloBatch(env,activityIntervals,batchIntervals,batchDecisions,loadMatrix[j]));

	}

	// End of new global constraint



	model.add(Cmax == IloMax(ends));

	// Upper Bound
	// model.add(Cmax <= 5000);

	// Provide initial solution? Assign job to its batch
	// Solve model adding constraint stating that a job is assigned to its batch, the optimize
	// Then use this solution as starting point in original model

	// Objective is to minimize makespan (maximum completion time across all jobs)
    IloObjective objective = IloMinimize(env,Cmax);
    model.add(objective);

	std::stringstream outputFile;

    outputFile << "C:/Users/Luigi/Desktop/Instances/output/batchCP_newOutput_" << x << "_instance_" << x << "_trial_" << trial << ".txt";
    std::string output = outputFile.str();

    //std::ifstream file("../../../examples/data/jobshop/jobshop_F2.txt");

    std::ofstream out(output);

	IloNum startClock = timer.getTime();

	// Flatten array of batch variables
	IloIntVarArray flatBatches(env);

	for (IloInt i = 0; i < nbJobs; i++) {

		for (IloInt j = 0; j < nbMachines; j++) {

			flatBatches.add(batchMatrix[i][j]);

		}

	}

	// Solve model
    IloCP cp(model);





	// Define one phase per job i.e. assign batchs of job 1 before moving on to job 2, and so on.
	// Maybe this would help the model find a solution more easily?
	// But could make it worse when trying to improve initial solution


    //cp.setParameter(IloCP::FailLimit, failLimit);
	cp.setParameter(IloCP::TimeMode, IloCP::ElapsedTime);
	cp.setParameter(IloCP::TimeLimit, 300);
	//cp.setParameter(IloCP::SearchType, IloCP::DepthFirst);
    cp.out() << "Instance \t: " << filename << std::endl;
    if (cp.solve(IloSearchPhase(env, flatBatches))) {
      cp.out() << "Makespan \t: " << cp.getObjValue() << std::endl;
	  out << "Makespan: " << cp.getObjValue() << ", Time: " << (cp.getTime());

	  for (IloInt i = 0; i < nbJobs; i++) {

		  for (IloInt j = 0; j < nbMachines; j++) {

			  cp.out() << " / " << cp.getValue(batchMatrix[i][j]);// << std::endl;

		  }

	  }


    } else {
      cp.out() << "No solution found."  << std::endl;
    }
  } catch(IloException& e){
    env.out() << " ERROR: " << e << std::endl;
  }

  }

  }

  }

  env.end();





  return 0;
}
