#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PhysicsListEM.hh"

#include "PrimaryGeneratorAction.hh"
#include "LHEPrimaryGeneratorAction.hh"
#include "SeededGeneratorAction.hh"

#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
#include "StackingAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc, char** argv) {
#ifdef G4VIS_USE
	std::cout << " -- G4VIS_USE is set " << std::endl;
#else
	std::cout << " -- G4VIS_USE is not set " << std::endl;
#endif

	// Choose the Random engine
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
	// User Verbose output class
	G4VSteppingVerbose::SetInstance(new SteppingVerbose);

	// Construct the default run manager
	G4RunManager * runManager = new G4RunManager;

	// Set mandatory initialization classes
	int version = DetectorConstruction::HGCAL_E26_TH;

	int model = DetectorConstruction::m_FULLSECTION;

	bool signal = true;
	std::string data = "";
	double steelThick= 0;
	if (argc > 2)
		version = atoi(argv[2]);
	if (argc > 3)
		model = atoi(argv[3]);
	if (argc > 4)
		signal = atoi(argv[4]);
	int status_x,status_y,seeds_x,seeds_y;
	if (argc > 7){
		status_x = atoi(argv[5]);
		status_y = atoi(argv[6]);
		seeds_x = atoi(argv[7]);
		seeds_y = atoi(argv[8]);

		ofstream outputFile("temp.rndm");
		Int_t line1 = status_x, line2 = status_y,
				line3 = seeds_x, line4 = seeds_y;
		outputFile << "Uvec\n" <<  line1 << "\n" << line2 << "\n"
				<< line3 << "\n" << line4;
	}

	CLHEP::HepRandom::restoreEngineStatus ("temp.rndm");

	std::cout << "-- Running version " << version << " model " << model
			<< std::endl;





	runManager->SetUserInitialization(
			new DetectorConstruction(version, model, steelThick));
	runManager->SetUserInitialization(new PhysicsList);

	// Set user action classes
	runManager->SetUserAction(new EventAction(signal));
	runManager->SetUserAction(new SteppingAction(data));
	runManager->SetUserAction(new StackingAction(data));

        if (data!="") {
            runManager->SetUserAction(new SeededGeneratorAction(model, data));
        }
        else {
            runManager->SetUserAction(new PrimaryGeneratorAction(model, signal, data));
        }
        runManager->Initialize();

	// Initialize visualization
#ifdef G4VIS_USE
	G4VisManager* visManager = new G4VisExecutive;
	visManager->Initialize();

#endif

	// Get the pointer to the User Interface manager
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	G4String fileName;
	if (argc > 1)
		fileName = argv[1];
	if (argc != 1)   // batch mode
			{
		std::cout << " ====================================== " << std::endl
				<< " ========  Running batch mode ========= " << std::endl
				<< " ====================================== " << std::endl;
		G4String command = "/control/execute ";
		UImanager->ApplyCommand(command + fileName);
	} else {
		std::cout << " ====================================== " << std::endl
				<< " ====  Running interactive display ==== " << std::endl
				<< " ====================================== " << std::endl;
#ifdef G4UI_USE
		G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
		UImanager->ApplyCommand("/control/execute vis.mac");
#endif
		if (ui->IsGUI())
		UImanager->ApplyCommand("/control/execute gui.mac");
		ui->SessionStart();
		delete ui;
#endif
	}

#ifdef G4VIS_USE
	delete visManager;
#endif
	delete runManager;

	return 0;
}
