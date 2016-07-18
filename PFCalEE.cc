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
	CLHEP::HepRandom::saveEngineStatus ("currentEvent.rndm");
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
	if (argc > 5)
		steelThick = std::stof(argv[5]);
	if (argc > 6)
		data = argv[6];
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
