#include "mainApp.h"
#include <wx/msgdlg.h> 
#include <wx/textdlg.h> 

wxIMPLEMENT_APP(RaytracerApp);



bool RaytracerApp::OnInit()
{
	if (!wxApp::OnInit())
		return false;



	// create the main frame window
	RaytracerFrame* frame = new RaytracerFrame();
	frame->SetSize(1200, 1000);

	renderPanel = new RenderPanel(frame, this);
	renderPanel->displayW = 600;
	renderPanel->displayH = 800;

	wxBoxSizer * statusbar_sizer = new wxBoxSizer(wxHORIZONTAL);
	progressBar = new wxGauge(frame->GetStatusBar(), wxID_ANY, 1000, wxDefaultPosition, wxDefaultSize, wxGA_HORIZONTAL);
	infoModel = new wxStaticText(frame->GetStatusBar(), wxID_ANY, "triangles: 0, vertices: 0, bvh leaves size: 0, mouse distance: 0.000");
	infoPerf = new wxStaticText(frame->GetStatusBar(), wxID_ANY, "Time per ray: 00.000 s                                           ");
	statusbar_sizer->Add(progressBar);
	statusbar_sizer->Add(infoPerf);
	statusbar_sizer->Add(infoModel);
	frame->GetStatusBar()->SetSizer(statusbar_sizer);

	m_bookCtrl = new wxNotebook(frame, wxID_ANY, wxDefaultPosition, wxSize(200, 100), 0);

	wxPanel *panelObject = new wxPanel(m_bookCtrl);
	wxBoxSizer * panelObject_sizer = new wxBoxSizer(wxVERTICAL);
	objectName = new wxTextCtrl(panelObject, 1000, "", wxDefaultPosition, wxDefaultSize);
	show_edges = new wxCheckBox(panelObject, EDGES_CHECKBOX, "Show Edges", wxDefaultPosition, wxDefaultSize);
	Connect(EDGES_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	interp_normals = new wxCheckBox(panelObject, INTERP_CHECKBOX, "Interp. Normals per Vtx", wxDefaultPosition, wxDefaultSize);
	Connect(INTERP_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);

	wxBoxSizer * transp_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* transp_text = new wxStaticText(panelObject, 9999, "Refraction index: ");
	transparent = new wxCheckBox(panelObject, TRANSPARENT_CHECKBOX, "Transparent", wxDefaultPosition, wxDefaultSize);
	Connect(TRANSPARENT_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	refractionIndex = new wxSpinCtrlDouble(panelObject, REFRACTION_INDEX, "1.5", wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 0.1, 10, 1.5, 0.05);
	Connect(REFRACTION_INDEX, wxEVT_SPINCTRLDOUBLE, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	transp_sizer->Add(transparent);
	transp_sizer->Add(transp_text);
	transp_sizer->Add(refractionIndex);

	flipnormals = new wxCheckBox(panelObject, FLIPNORMALS_CHECKBOX, "Flip normals for transparency", wxDefaultPosition, wxDefaultSize);
	Connect(FLIPNORMALS_CHECKBOX, wxEVT_COMMAND_CHECKBOX_CLICKED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);


	m_AlbedoFile = new wxListCtrl(panelObject, ALBEDO_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol;
	itCol.SetId(0);
	itCol.SetText("Albedo Texture");
	itCol.SetWidth(200);
	wxListItem itCol2;
	itCol2.SetId(1);
	itCol2.SetText("Albedo Color");
	itCol2.SetWidth(200);
	m_AlbedoFile->InsertColumn(0, itCol);
	m_AlbedoFile->InsertColumn(1, itCol2);
	m_AlbedoFile->SetDropTarget(new DnDAlbedoFile(m_AlbedoFile, frame));
	Connect(ALBEDO_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClick), NULL, frame);
	Connect(ALBEDO_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClick), NULL, frame);

	m_SpecularFile = new wxListCtrl(panelObject, SPECULAR_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol4;
	itCol4.SetId(0);
	itCol4.SetText("Specular Map");
	itCol4.SetWidth(200);
	itCol2.SetId(1);
	itCol2.SetText("Specular Color");
	itCol2.SetWidth(200);
	m_SpecularFile->InsertColumn(0, itCol4);
	m_SpecularFile->InsertColumn(1, itCol2);
	m_SpecularFile->SetDropTarget(new DnDSpecularFile(m_SpecularFile, frame));
	Connect(SPECULAR_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickSpecular), NULL, frame);
	Connect(SPECULAR_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickSpecular), NULL, frame);

	m_RoughnessFile = new wxListCtrl(panelObject, ROUGHNESS_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol5;
	itCol5.SetId(0);
	itCol5.SetText("Phong Exponent Map");
	itCol5.SetWidth(200);
	itCol2.SetId(1);
	itCol2.SetText("Exponent Values");
	itCol2.SetWidth(200);
	m_RoughnessFile->InsertColumn(0, itCol5);
	m_RoughnessFile->InsertColumn(1, itCol2);
	m_RoughnessFile->SetDropTarget(new DnDRoughnessFile(m_RoughnessFile, frame));
	Connect(ROUGHNESS_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickRoughness), NULL, frame);
	Connect(ROUGHNESS_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickRoughness), NULL, frame);

	m_NormalFile = new wxListCtrl(panelObject, NORMAL_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	itCol2.SetId(0);
	itCol2.SetText("Normal Map");
	itCol2.SetWidth(300);
	m_NormalFile->InsertColumn(0, itCol2);
	m_NormalFile->SetDropTarget(new DnDNormalFile(m_NormalFile, frame));
	Connect(NORMAL_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickNormal), NULL, frame);
	Connect(NORMAL_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickNormal), NULL, frame);

	//albedoColorPicker = new wxColourPickerCtrl(panelObject, ALBEDO_COLORPICKER, wxColour(255,255,255), wxDefaultPosition, wxDefaultSize, wxCLRP_USE_TEXTCTRL | wxCLRP_SHOW_LABEL);
	//Connect(ALBEDO_COLORPICKER, wxEVT_COLOURPICKER_CHANGED, wxCommandEventHandler(RenderPanel<double>::update_parameters_and_render), NULL, renderPanel);

	m_AlphaFile = new wxListCtrl(panelObject, ALPHA_FILES, wxDefaultPosition, wxSize(-1, 100), wxLC_REPORT);
	wxListItem itCol3;
	itCol3.SetId(0);
	itCol3.SetText("Alpha Map");
	itCol3.SetWidth(300);
	m_AlphaFile->InsertColumn(0, itCol3);
	m_AlphaFile->SetDropTarget(new DnDAlphaFile(m_AlphaFile, frame));
	Connect(ALPHA_FILES, wxEVT_LIST_ITEM_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickAlpha), NULL, frame);
	Connect(ALPHA_FILES, wxEVT_LIST_COL_RIGHT_CLICK, wxListEventHandler(RaytracerFrame::OnListRightClickAlpha), NULL, frame);

	/*wxBoxSizer * ks_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* ks_text = new wxStaticText(panelObject, 9999, "Ks : ");
	ks_slider = new wxSlider(panelObject, KS_SLIDER, 0, 0, 100, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL);
	Connect(KS_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel<double>::update_parameters_and_render), NULL, renderPanel);
	ks_sizer->Add(ks_text, 0, wxEXPAND);
	ks_sizer->Add(ks_slider, 1, wxEXPAND);*/

	deleteObject = new wxButton(panelObject, DELETE_OBJECT, "Delete", wxDefaultPosition, wxDefaultSize);
	Connect(DELETE_OBJECT, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::delete_object), NULL, renderPanel);

	panelObject_sizer->Add(objectName, 0, wxEXPAND);
	panelObject_sizer->Add(show_edges, 0, wxEXPAND);
	panelObject_sizer->Add(interp_normals, 0, wxEXPAND);
	panelObject_sizer->Add(transp_sizer, 0, wxEXPAND);
	panelObject_sizer->Add(flipnormals, 0, wxEXPAND);
	panelObject_sizer->Add(m_AlbedoFile, 0, wxEXPAND);
	//panelObject_sizer->Add(albedoColorPicker, 0, wxEXPAND);
	panelObject_sizer->Add(m_SpecularFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_RoughnessFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_NormalFile, 0, wxEXPAND);
	panelObject_sizer->Add(m_AlphaFile, 0, wxEXPAND);
	//panelObject_sizer->Add(ks_sizer, 0, wxEXPAND);
	panelObject_sizer->Add(deleteObject, 0, wxEXPAND);
	
	panelObject->SetSizer(panelObject_sizer);


	m_bookCtrl->AddPage(panelObject, wxT("Object"), false);


	wxPanel *panelRenderer = new wxPanel(m_bookCtrl);
	wxBoxSizer * panelRenderer_sizer = new wxBoxSizer(wxVERTICAL);

	wxBoxSizer * rendersize_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* renderwidth_text = new wxStaticText(panelRenderer, 9999 - 1, "W: ");
	renderwidth = new wxSpinCtrl(panelRenderer, RENDER_WIDTH, wxString("1000"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 2, 10000, 1000);	
	Connect(RENDER_WIDTH, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	wxStaticText* renderheight_text = new wxStaticText(panelRenderer, 9999 - 1, "H: ");
	renderheight = new wxSpinCtrl(panelRenderer, RENDER_HEIGHT, wxString("800"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 2, 10000, 800);
	Connect(RENDER_HEIGHT, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	rendersize_sizer->Add(renderwidth_text, 0, wxEXPAND);
	rendersize_sizer->Add(renderwidth, 1, wxEXPAND);
	rendersize_sizer->Add(renderheight_text, 0, wxEXPAND);
	rendersize_sizer->Add(renderheight, 1, wxEXPAND);

	wxBoxSizer * nbrays_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* nbrays_text = new wxStaticText(panelRenderer, 9999 - 1, "nb paths per pixel: ");
	nbrays = new wxSpinCtrl(panelRenderer, NBRAYS, wxString("1000"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 10000, 1000);
	Connect(NBRAYS, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	nbrays_sizer->Add(nbrays_text, 0, wxEXPAND);
	nbrays_sizer->Add(nbrays, 1, wxEXPAND);

	wxBoxSizer * fov_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fov_text = new wxStaticText(panelRenderer, 9999, "fov (deg): ");
	fov_slider = new wxSlider(panelRenderer, FOV_SLIDER, 35, 2, 179, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL|wxSL_LABELS);
	Connect(FOV_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fov_sizer->Add(fov_text, 0, wxEXPAND);
	fov_sizer->Add(fov_slider, 1, wxEXPAND);

	wxBoxSizer * aperture_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* aperture_text = new wxStaticText(panelRenderer, 9999-1, "aperture (mm): ");
	aperture_slider = new wxSlider(panelRenderer, APERTURE_SLIDER, 10, 0, 3000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(APERTURE_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	aperture_sizer->Add(aperture_text, 0, wxEXPAND);
	aperture_sizer->Add(aperture_slider, 1, wxEXPAND);

	wxBoxSizer * focus_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* focus_text = new wxStaticText(panelRenderer, 9999 - 1, "focus distance (cm): ");
	focus_slider = new wxSlider(panelRenderer, FOCUS_SLIDER, 5000, 0, 10000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(FOCUS_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	focus_sizer->Add(focus_text, 0, wxEXPAND);
	focus_sizer->Add(focus_slider, 1, wxEXPAND);

	wxBoxSizer * filter_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* filter_text = new wxStaticText(panelRenderer, 9999 - 1, "filter width: ");
	filter_slider = new wxSlider(panelRenderer, FILTER_SLIDER, 5, 1, 100, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(FILTER_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	filter_sizer->Add(filter_text, 0, wxEXPAND);
	filter_sizer->Add(filter_slider, 1, wxEXPAND);


	wxBoxSizer * envmap_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* envmap_text = new wxStaticText(panelRenderer, 9999 - 1, "envmap: ");
	envmapName = new wxTextCtrl(panelRenderer, 1000, "", wxDefaultPosition, wxDefaultSize);
	envmapName->SetDropTarget(new DnDEnvmapFile(envmapName, frame));
	envmap_sizer->Add(envmap_text, 0, wxEXPAND);
	envmap_sizer->Add(envmapName, 1, wxEXPAND);

	wxBoxSizer * envmapintensity_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* envmapintensity_text = new wxStaticText(panelRenderer, 9999 - 1, "envmap intensity: ");
	envmapintensity_slider = new wxSlider(panelRenderer, ENVMAPINTENSITY_SLIDER, 100, 0, 1000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(ENVMAPINTENSITY_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	envmapintensity_sizer->Add(envmapintensity_text, 0, wxEXPAND);
	envmapintensity_sizer->Add(envmapintensity_slider, 1, wxEXPAND);

	wxBoxSizer * lightintensity_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* lightintensity_text = new wxStaticText(panelRenderer, 9999 - 1, "light intensity: ");
	lightintensity_slider = new wxSlider(panelRenderer, LIGHTINTENSITY_SLIDER, 100, 0, 1000, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(LIGHTINTENSITY_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	lightintensity_sizer->Add(lightintensity_text, 0, wxEXPAND);
	lightintensity_sizer->Add(lightintensity_slider, 1, wxEXPAND);

	wxBoxSizer * bounces_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* bounces_text = new wxStaticText(panelRenderer, 9999 - 1, "light bounces: ");
	bounces = new wxSpinCtrl(panelRenderer, BOUNCES_SPIN, wxString("3"), wxDefaultPosition, wxDefaultSize, wxSP_ARROW_KEYS, 1, 10, 3);
	Connect(BOUNCES_SPIN, wxEVT_SPINCTRL, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	bounces_sizer->Add(bounces_text, 0, wxEXPAND);
	bounces_sizer->Add(bounces, 1, wxEXPAND);

	launchRender = new wxButton(panelRenderer, LAUNCH_RENDER, "Offline render (freezes app)", wxDefaultPosition, wxDefaultSize);
	Connect(LAUNCH_RENDER, wxEVT_BUTTON, wxCommandEventHandler(RenderPanel::launch_render), NULL, renderPanel);
	
	panelRenderer_sizer->Add(rendersize_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(nbrays_sizer, 0, wxEXPAND);	
	panelRenderer_sizer->Add(fov_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(aperture_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(focus_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(filter_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(envmap_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(envmapintensity_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(lightintensity_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(bounces_sizer, 0, wxEXPAND);
	panelRenderer_sizer->Add(launchRender, 0, wxEXPAND);
	panelRenderer->SetSizer(panelRenderer_sizer);

	m_bookCtrl->AddPage(panelRenderer, wxT("Renderer"), false);



	wxPanel *panelFog = new wxPanel(m_bookCtrl);
	wxBoxSizer * panelFog_sizer = new wxBoxSizer(wxVERTICAL);

	wxBoxSizer * fogdensity_sizer = new wxBoxSizer(wxHORIZONTAL);
	wxStaticText* fogdensity_text = new wxStaticText(panelFog, 9999 - 1, "fog density : ");
	fogdensity_slider = new wxSlider(panelFog, FOGDENSITY_SLIDER, 0, 0, 100, wxDefaultPosition, wxDefaultSize, wxSL_HORIZONTAL | wxSL_LABELS);
	Connect(FOGDENSITY_SLIDER, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	fogdensity_sizer->Add(fogdensity_text, 0, wxEXPAND);
	fogdensity_sizer->Add(fogdensity_slider, 1, wxEXPAND);

	uniformFogRadio = new wxRadioButton(panelFog, UNIFORMFOG_RADIO, "Uniform fog", wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
	expFogRadio = new wxRadioButton(panelFog, EXPFOG_RADIO, "Exponential fog", wxDefaultPosition, wxDefaultSize, 0);
	Connect(UNIFORMFOG_RADIO, wxEVT_RADIOBUTTON, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);
	Connect(UNIFORMFOG_RADIO, wxEVT_RADIOBUTTON, wxCommandEventHandler(RenderPanel::update_parameters_and_render), NULL, renderPanel);

	panelFog_sizer->Add(fogdensity_sizer, 0, wxEXPAND);
	panelFog_sizer->Add(uniformFogRadio, 0, wxEXPAND);
	panelFog_sizer->Add(expFogRadio, 0, wxEXPAND);
	panelFog->SetSizer(panelFog_sizer);

	m_bookCtrl->AddPage(panelFog, wxT("Fog"), false);



	m_bookCtrl->SetSelection(0);


	wxBoxSizer * sizer = new wxBoxSizer(wxHORIZONTAL);
	frame->SetSizer(sizer);
	sizer->Add(renderPanel, 2, wxEXPAND);
	sizer->Add(m_bookCtrl, 1, wxEXPAND);


	colPicker = new wxColourDialog(frame);
	texOpenDlg = new wxFileDialog(frame, _("Open texture file"), "", "", "All files (*.*)|*.*", wxFD_OPEN | wxFD_FILE_MUST_EXIST);

	render_loop_on = false;
	frame->Show();
	activateRenderLoop(true);

    return true;
}


void RaytracerApp::activateRenderLoop(bool on)
{
	renderPanel->SetDoubleBuffered(true);
	if (on && !render_loop_on)
	{
		Connect(wxID_ANY, wxEVT_IDLE, wxIdleEventHandler(RaytracerApp::OnIdle));
		Connect(wxID_ANY, wxEVT_ERASE_BACKGROUND, wxEraseEventHandler(RaytracerApp::OnEraseBackGround));
		render_loop_on = true;
	} else if (!on && render_loop_on) {
		Disconnect(wxEVT_IDLE, wxIdleEventHandler(RaytracerApp::OnIdle));
		render_loop_on = false;
	}
}

RaytracerFrame::RaytracerFrame()
        : wxFrame(NULL, wxID_ANY, wxT("Pathtracer by N.Bonneel"),
                  wxPoint(10, 100))

{
    // frame icon and status bar
    SetIcon(wxICON(sample));

    CreateStatusBar();

    // construct menu
    wxMenu *file_menu = new wxMenu();
    file_menu->Append(ID_SAVE, wxT("&Save scene as..."));
	Connect(ID_SAVE, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::SaveAs), NULL, this);
	file_menu->Append(ID_OPEN, wxT("&Open scene..."));
	Connect(ID_OPEN, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::Open), NULL, this);
	file_menu->Append(ID_SAVEIMAGE, wxT("&Save current image..."));
	Connect(ID_SAVEIMAGE, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::SaveImage), NULL, this);
	file_menu->Append(ID_EXPORT_MTL, wxT("&Export .mtl (selected object)..."));
	Connect(ID_EXPORT_MTL, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::ExportMtl), NULL, this);

	wxMenu *info_menu = new wxMenu();
	info_menu->Append(ID_MESHINFO, wxT("&Mesh information..."));
	Connect(ID_MESHINFO, wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::ShowMeshInfo), NULL, this);

	wxMenuBar* menu_bar = new wxMenuBar();
	menu_bar->Append(file_menu, wxT("&File"));	
	menu_bar->Append(info_menu, wxT("&Info"));

    SetMenuBar(menu_bar);

    SetBackgroundColour(*wxWHITE); // labels read better on this background
	
    this->SetDropTarget(new DnDFile(this));	

	this->SetSize(1024, 1024);
	Centre();

}




void RenderPanel::update_parameters_and_render(wxCommandEvent& event) {

	if (selected_object < 0) return;
	if (selected_object >= raytracer.s.objects.size()) return;

	stop_render();

	raytracer.W = raytracer_app->renderwidth->GetValue();
	raytracer.H = raytracer_app->renderheight->GetValue();
	raytracer.nrays = raytracer_app->nbrays->GetValue();
	raytracer.s.objects[selected_object]->display_edges = raytracer_app->show_edges->IsChecked();
	raytracer.s.objects[selected_object]->interp_normals = raytracer_app->interp_normals->IsChecked();
	raytracer.s.objects[selected_object]->transparent = raytracer_app->transparent->IsChecked();
	raytracer.cam.fov = raytracer_app->fov_slider->GetValue() * M_PI / 180.;
	raytracer.cam.aperture = raytracer_app->aperture_slider->GetValue() / 1000.;
	raytracer.cam.focus_distance = raytracer_app->focus_slider->GetValue() / 100.;
	raytracer.sigma_filter = raytracer_app->filter_slider->GetValue() / 10.;
	//wxColour alb = raytracer_app->albedoColorPicker->GetColour();
	//raytracer.s.objects[selected_object]->albedo = Vector(alb.Red()/255., alb.Green()/255., alb.Blue()/255.);

	//raytracer.s.objects[selected_object]->ks = raytracer_app->ks_slider->GetValue() / 100.;
	raytracer.s.fog_density = raytracer_app->fogdensity_slider->GetValue() / 100.;

	raytracer.s.intensite_lumiere = raytracer_app->lightintensity_slider->GetValue() / 100. * 1000000000 * 4.*M_PI / (4.*M_PI*raytracer.s.lumiere->R*raytracer.s.lumiere->R*M_PI);
	raytracer.s.envmap_intensity = raytracer_app->envmapintensity_slider->GetValue() / 100.;

	raytracer.s.fog_type = raytracer_app->uniformFogRadio->GetValue() ? 0 : 1;
	raytracer.nb_bounces = raytracer_app->bounces->GetValue();

	raytracer.s.objects[selected_object]->flip_normals = raytracer_app->flipnormals->IsChecked();
	raytracer.s.objects[selected_object]->refr_index = raytracer_app->refractionIndex->GetValue();



	start_render();
}
void RenderPanel::update_textures_and_render(wxCommandEvent& event) {

	if (selected_object < 0) return;
	if (selected_object >= raytracer.s.objects.size()) return;

	stop_render();

	raytracer.s.objects[selected_object]->display_edges = raytracer_app->show_edges->IsChecked();

	start_render();
}

void RenderPanel::update_gui() {

	if ((selected_object < 0) || (selected_object >= raytracer.s.objects.size())) {
		raytracer_app->objectName->SetLabelText(" ");
		return;
	}

	Geometry* g = dynamic_cast<Geometry*>(raytracer.s.objects[selected_object]);
	if (g) {
		std::ostringstream os;
		os << "triangles: " << g->indices.size() << ", vertices: " << g->vertices.size() << ", bvh leaves size: " << g->max_bvh_triangles << ", bvh max depth: " << g->bvh_depth << ", bvh avg depth:" << g->bvh_avg_depth << ", bvh nb nodes: " << g->bvh_nb_nodes << ", mouse distance: " << selected_object_t;
		raytracer_app->infoModel->SetLabelText(os.str().c_str());
	} else {
		std::ostringstream os;
		os << "mouse distance: " << selected_object_t;
		raytracer_app->infoModel->SetLabelText(os.str().c_str());
	}

	raytracer_app->renderwidth->SetValue(raytracer.W);
	raytracer_app->renderheight->SetValue(raytracer.H);
	raytracer_app->nbrays->SetValue(raytracer.nrays);
	raytracer_app->objectName->SetLabelText(raytracer.s.objects[selected_object]->name);
	raytracer_app->show_edges->SetValue(raytracer.s.objects[selected_object]->display_edges);
	raytracer_app->interp_normals->SetValue(raytracer.s.objects[selected_object]->interp_normals);
	raytracer_app->fov_slider->SetValue(raytracer.cam.fov * 180. / M_PI);
	raytracer_app->aperture_slider->SetValue(raytracer.cam.aperture * 1000);
	raytracer_app->focus_slider->SetValue(raytracer.cam.focus_distance*100.);
	//raytracer_app->ks_slider->SetValue(raytracer.s.objects[selected_object]->ks*100.);
	raytracer_app->filter_slider->SetValue(raytracer.sigma_filter*10.);
	raytracer_app->bounces->SetValue(raytracer.nb_bounces);
	raytracer_app->transparent->SetValue(raytracer.s.objects[selected_object]->transparent);
	raytracer_app->flipnormals->SetValue(raytracer.s.objects[selected_object]->flip_normals);
	raytracer_app->refractionIndex->SetValue(raytracer.s.objects[selected_object]->refr_index);

	//Vector alb = raytracer.s.objects[selected_object]->albedo;
	//raytracer_app->albedoColorPicker->SetColour(wxColour(alb[0] * 255., alb[1] * 255., alb[2] * 255.));

	raytracer_app->m_AlbedoFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->textures.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->textures[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->textures[i].filename;
		long index = raytracer_app->m_AlbedoFile->InsertItem(i, item);
		raytracer_app->m_AlbedoFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_AlbedoFile->SetItem(index, 1, raytracer.s.objects[selected_object]->textures[i].multiplier.toColorStr(), -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_AlbedoFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	raytracer_app->m_NormalFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->normal_map.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->normal_map[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->normal_map[i].filename;
		long index = raytracer_app->m_NormalFile->InsertItem(i, item);
		raytracer_app->m_NormalFile->SetItem(index, 0, txt, -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_NormalFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}
	raytracer_app->envmapName->SetLabelText(((Sphere*)raytracer.s.objects[1])->envmapfilename);

	raytracer_app->m_AlphaFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->alphamap.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->alphamap[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->alphamap[i].filename;
		long index = raytracer_app->m_AlphaFile->InsertItem(i, item);
		raytracer_app->m_AlphaFile->SetItem(index, 0, txt, -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_AlphaFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	raytracer_app->m_SpecularFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->specularmap.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->specularmap[i].filename);
		wxListItem item;
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->specularmap[i].filename;
		long index = raytracer_app->m_SpecularFile->InsertItem(i, item);
		raytracer_app->m_SpecularFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_SpecularFile->SetItem(index, 1, raytracer.s.objects[selected_object]->specularmap[i].multiplier.toColorStr(), -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_SpecularFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	raytracer_app->m_RoughnessFile->DeleteAllItems();
	for (int i = 0; i < raytracer.s.objects[selected_object]->roughnessmap.size(); i++) {
		std::string txt = extractFileName(raytracer.s.objects[selected_object]->roughnessmap[i].filename);
		wxListItem item;
		item.SetText(txt);
		item.SetId(i);
		std::string filename = raytracer.s.objects[selected_object]->roughnessmap[i].filename;
		long index = raytracer_app->m_RoughnessFile->InsertItem(i, item);
		raytracer_app->m_RoughnessFile->SetItem(index, 0, txt, -1);
		raytracer_app->m_RoughnessFile->SetItem(index, 1, raytracer.s.objects[selected_object]->roughnessmap[i].multiplier.toColorStr(), -1);
		if (filename[0] != 'N' && filename[1] != 'u' && !file_exists(filename.c_str()))
			raytracer_app->m_RoughnessFile->SetItemBackgroundColour(index, wxColour(255, 0, 0));
	}

	if (g) {
		int id = g->indices[selected_tri].group;
		raytracer_app->m_AlbedoFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_NormalFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_AlphaFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_SpecularFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);
		raytracer_app->m_RoughnessFile->SetItemState(id, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED, wxLIST_STATE_SELECTED | wxLIST_STATE_FOCUSED);

		raytracer_app->m_AlbedoFile->EnsureVisible(id);
		raytracer_app->m_NormalFile->EnsureVisible(id);
		raytracer_app->m_AlphaFile->EnsureVisible(id);
		raytracer_app->m_SpecularFile->EnsureVisible(id);
		raytracer_app->m_RoughnessFile->EnsureVisible(id);
	}

	raytracer_app->fogdensity_slider->SetValue(raytracer.s.fog_density*100.);
	raytracer_app->uniformFogRadio->SetValue(raytracer.s.fog_type == 0);

	double factor = 1. / 100. * 1000000000 * 4.*M_PI / (4.*M_PI*raytracer.s.lumiere->R*raytracer.s.lumiere->R*M_PI);
	raytracer_app->lightintensity_slider->SetValue(raytracer.s.intensite_lumiere / factor);
	raytracer_app->envmapintensity_slider->SetValue(raytracer.s.envmap_intensity * 100);
}

void RenderPanel::render(wxDC& dc)
{
	if (cur_img.size() == 0) return;
	if (raytracer.stopped) return;
	raytracer_app->progressBar->SetValue(raytracer.current_nb_rays / (double)raytracer.nrays*1000.);
	std::ostringstream os;
	os << "Time per ray: " << raytracer.curTimePerFrame / 1000. << " s";
	raytracer_app->infoPerf->SetLabelText(os.str().c_str());

	gamma_corrected_image.resize(raytracer.W*raytracer.H * 3);
	extrapolated_image = raytracer.imagedouble;
	/*computed = raytracer.computed;
	computed2 = raytracer.computed;
	double lastR = 0, lastG = 0, lastB = 0;
	int offsetsi[4] = { -1, 0, 0, 1 };
	int offsetsj[4] = { 0, -1, 1, 0 };
	int nbchanged;
	for (int iter = 0; iter < 10; iter++) {   // Voronoi : lent et moche
	nbchanged = 0;
	#pragma omp parallel for
	for (int i = 0; i < raytracer.H; i++) {
	for (int j = 0; j < raytracer.W; j++) {
	if (!computed2[i*raytracer.W + j]) {
	//for (int k = std::max(0, i - 1); k <= std::min(H - 1, i + 1); k++) {
	//	for (int l = std::max(0, j - 1); l <= std::min(W - 1, j + 1); l++) {
	int rand_shuff = rand();
	for (int id=0; id<4; id++) {
	int k = i + offsetsi[(id+ rand_shuff)%4]; if (k < 0 || k>= raytracer.H) continue;
	int l = j + offsetsj[(id + rand_shuff) % 4]; if (l < 0 || l >= raytracer.W) continue;
	if (computed2[k*raytracer.W + l]) {
	extrapolated_image[(i*raytracer.W + j) * 3 + 0] = extrapolated_image[(k*raytracer.W + l) * 3 + 0];
	extrapolated_image[(i*raytracer.W + j) * 3 + 1] = extrapolated_image[(k*raytracer.W + l) * 3 + 1];
	extrapolated_image[(i*raytracer.W + j) * 3 + 2] = extrapolated_image[(k*raytracer.W + l) * 3 + 2];
	computed[i*raytracer.W + j] = true;
	nbchanged++;
	break;
	}
	}
	}
	}
	}
	if (nbchanged == 0) break;
	computed2 = computed;
	}*/

	for (int i = 0; i < raytracer.H; i++) {
		for (int j = 0; j < raytracer.W; j++) {
			//if (!raytracer.computed[i*raytracer.W + j]) {
			//if (j< raytracer.W/2) {
			if (raytracer.sample_count[i*raytracer.W + j] <= 5) {
				double alpha = raytracer.sample_count[i*raytracer.W + j] / 6.;

				int ix = j / 16;
				double fx = (j / 16.) - ix;

				int iy = i / 16;
				double fy = (i / 16.) - iy;

				for (int k = 0; k < 3; k++) {
					//double lowval = raytracer.imagedouble_lowres[(iy * raytracer.Wlr + ix) * 3 + k];  // nearest neighbor interp.
					double lowval;
					if (ix < raytracer.Wlr - 1 && iy < raytracer.Hlr - 1)
						lowval = (raytracer.imagedouble_lowres[(iy*raytracer.Wlr + ix) * 3 + k] * (1 - fx) + raytracer.imagedouble_lowres[(iy*raytracer.Wlr + ix + 1) * 3 + k] * fx) * (1 - fy)
						+ (raytracer.imagedouble_lowres[((iy + 1)*raytracer.Wlr + ix) * 3 + k] * (1 - fx) + raytracer.imagedouble_lowres[((iy + 1)*raytracer.Wlr + ix + 1) * 3 + k] * fx) * fy;
					else
						lowval = raytracer.imagedouble_lowres[(iy*raytracer.Wlr + ix) * 3 + k];
					extrapolated_image[(i*raytracer.W + j) * 3 + k] = extrapolated_image[(i*raytracer.W + j) * 3 + k] * alpha + (1. - alpha)*lowval;
				}
			}
		}
	}
	//os << "Nbchanged: " << nbchanged;
	//raytracer_app->infoPerf->SetLabelText(os.str().c_str());

#pragma omp parallel for 
	for (int i = 0; i < raytracer.H; i++) {
		for (int j = 0; j < raytracer.W; j++) {
			double normalization = 1. / (raytracer.sample_count[i*raytracer.W + j] + 1.);
			/*gamma_corrected_image[(i*raytracer.W + j) * 3 + 0] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 0] / (raytracer.current_nb_rays + 1)), 1 / 2.2));   // rouge
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 1] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 1] / (raytracer.current_nb_rays + 1)), 1 / 2.2)); // vert
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 2] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 2] / (raytracer.current_nb_rays + 1)), 1 / 2.2)); // bleu*/
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 0] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 0] * normalization), 1 / 2.2));   // rouge
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 1] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 1] * normalization), 1 / 2.2)); // vert
			gamma_corrected_image[(i*raytracer.W + j) * 3 + 2] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 2] * normalization), 1 / 2.2)); // bleu

		}
	}

	static int nbcalls = -1; nbcalls++;
	if (nbcalls == 0 || screenImage.GetWidth() != raytracer.W || screenImage.GetHeight() != raytracer.H) {
		screenImage = wxImage(raytracer.W, raytracer.H, &(gamma_corrected_image[0]), true);
	}

	bmpBuf = wxBitmap(screenImage, dc);


	double scale_x = (double)displayW / raytracer.W;
	double scale_y = (double)displayH / raytracer.H;
	dc.SetUserScale(scale_x, scale_y);
	dc.DrawBitmap(bmpBuf, 0, 0);
	dc.SetUserScale(1.0, 1.0);
}


void RaytracerFrame::SaveAs(wxCommandEvent &evt) {
	wxFileDialog
		saveFileDialog(this, _("Save SCN file"), "", "",
			"SCN files (*.scn)|*.scn", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (saveFileDialog.ShowModal() == wxID_CANCEL)
		return;     

	render_panel->raytracer.save_scene(saveFileDialog.GetPath());	
}


void RaytracerFrame::SaveImage(wxCommandEvent &evt) {
	wxFileDialog
		saveFileDialog(this, _("Save Image file"), "", "",
			"Image files (*.*)|*.*", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (saveFileDialog.ShowModal() == wxID_CANCEL)
		return;

	int nbr = render_panel->raytracer.current_nb_rays;
	int W = render_panel->raytracer.W;
	int H = render_panel->raytracer.H;
	std::vector<unsigned char> image(W*H * 3);
#pragma omp parallel for 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(render_panel->raytracer.imagedouble[((H - i - 1)*W + j) * 3 + 0] / (nbr + 1.), 1 / 2.2)));   // rouge
			image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(render_panel->raytracer.imagedouble[((H - i - 1)*W + j) * 3 + 1] / (nbr + 1.), 1 / 2.2))); // vert
			image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(render_panel->raytracer.imagedouble[((H - i - 1)*W + j) * 3 + 2] / (nbr + 1.), 1 / 2.2))); // bleu
		}
	}
	save_img(saveFileDialog.GetPath().c_str(), &image[0], W, H);
}


void RaytracerFrame::ExportMtl(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	Geometry* g = dynamic_cast<Geometry*>(render_panel->raytracer.s.objects[obj_id]);
	if (!g) return;

	wxFileDialog
		saveFileDialog(this, _("Export material"), "", "",
			"MTL files (*.mtl)|*.mtl", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
	if (saveFileDialog.ShowModal() == wxID_CANCEL)
		return;

	g->exportMTL(saveFileDialog.GetPath().c_str());
}

void RaytracerFrame::Open(wxCommandEvent &evt) {
	wxFileDialog
		openFileDialog(this, _("Open SCN file"), "", "",
			"SCN files (*.scn)|*.scn", wxFD_OPEN | wxFD_FILE_MUST_EXIST);
	if (openFileDialog.ShowModal() == wxID_CANCEL)
		return;

	render_panel->stop_render();
	render_panel->raytracer.load_scene(openFileDialog.GetPath());
	render_panel->update_gui();
	render_panel->start_render();
}


void RaytracerFrame::ShowMeshInfo(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	Geometry* g = dynamic_cast<Geometry*>(render_panel->raytracer.s.objects[obj_id]);
	if (!g) return;

	int nbE, nbNonM, nbBoundaries;
	int nbC = g->getNbConnected(nbE, nbNonM, nbBoundaries);
	int Euler = (int)g->vertices.size() - nbE + (int)g->indices.size();
	int nbRealEdges, nbOthers, nbTri;
	g->findQuads(nbTri, nbOthers, nbRealEdges);

	std::ostringstream os;
	os << "Object ID: " << obj_id << std::endl;
	os << "Nb Tessellated Triangles: " << g->indices.size() << std::endl;
	os << "Nb Triangles in Original mesh: " << nbTri << std::endl;
	os << "Nb Quads (or other poly) in Original mesh: " << nbOthers << std::endl;
	os << "Nb Vertices: " << g->vertices.size() << std::endl;
	os << "Nb Tessellated Triangle Edges: " << nbE << std::endl;
	os << "Nb Edges in Original mesh: " << nbRealEdges << std::endl;
	os << "Nb Connected Components: " << nbC << std::endl;
	os << "Non Manifold Edges: " << nbNonM << std::endl;
	os << "Boundary edges: " << nbBoundaries << std::endl;
	os << "Euler Number: " << Euler << std::endl;
	os << "Genus: " << (2 - Euler)/2 << std::endl;
	os << " BVH max depth: " << g->bvh_depth << std::endl;
	os << " BVH avg depth: " << g->bvh_avg_depth << std::endl;
	os << " BVH nb nodes: " << g->bvh_nb_nodes << std::endl;
	os << " BVH max triangles per leaf: " << g->max_bvh_triangles << std::endl;

	wxMessageBox(os.str().c_str(), "Mesh information");


}

void RaytracerFrame::OnPopupClick(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = (int)static_cast<wxMenu *>(evt.GetEventObject())->GetClientData();
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_AlbedoFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);

	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_AlbedoFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_texture(firstSel);
		}
		render_panel->start_render();
		break;	
	case ID_MOVEUP:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_textures(item_id, item_id-1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN:
		if (item_id >= render_panel->raytracer_app->m_AlbedoFile->GetItemCount()-1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_textures(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE:
		if (item_id >= render_panel->raytracer_app->m_AlbedoFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->textures[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_ADDWHITE:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_texture(Vector(col.Red()/255., col.Green() / 255., col.Blue() / 255.));
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_COLOR:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_texture(Vector(col.Red() / 255., col.Green() / 255., col.Blue() / 255.), item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();			
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_texture(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}
void RaytracerFrame::OnPopupClickSpecular(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = (int)static_cast<wxMenu *>(evt.GetEventObject())->GetClientData();
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_SpecularFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);

	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE_SPECULAR:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_SpecularFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_specular(firstSel);
		}		
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE_SPECULAR:
		if (item_id >= render_panel->raytracer_app->m_AlbedoFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->specularmap[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_MOVEUP_SPECULAR:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_specular(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_SPECULAR:
		if (item_id >= render_panel->raytracer_app->m_SpecularFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_specular(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_SPECULAR:
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_specular(Vector(col.Red() / 255., col.Green() / 255., col.Blue() / 255.));
			render_panel->start_render();			
		}
		break;
	case ID_CHANGE_COLOR_SPECULAR:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_specular(Vector(col.Red() / 255., col.Green() / 255., col.Blue() / 255.), item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_SPECULAR:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_specularmap(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}
void RaytracerFrame::OnPopupClickRoughness(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = (int)static_cast<wxMenu *>(evt.GetEventObject())->GetClientData();
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_RoughnessFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);

	switch (evt.GetId()) {
	case ID_DELETE_ROUGHNESS:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_RoughnessFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_roughness(firstSel);
		}		
		render_panel->start_render();
		break;
	case ID_REMOVE_TEXTURE_ROUGHNESS:
		if (item_id >= render_panel->raytracer_app->m_AlbedoFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->roughnessmap[item_id].clear_texture();
		render_panel->start_render();
		break;
	case ID_MOVEUP_ROUGHNESS:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_roughness(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_ROUGHNESS:
		if (item_id >= render_panel->raytracer_app->m_RoughnessFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_roughness(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_ROUGHNESS:
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_roughness(Vector(col.Red(), col.Green(), col.Blue()));
			render_panel->start_render();			
		}
		break;
	case ID_CHANGE_COLOR_ROUGHNESS:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_roughness(Vector(col.Red(), col.Green(), col.Blue()), item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_ROUGHNESS:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_roughnessmap(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}
void RaytracerFrame::OnPopupClickNormal(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = (int)static_cast<wxMenu *>(evt.GetEventObject())->GetClientData();

	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_NormalFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE_NORMAL:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_NormalFile->GetNextItem(itemIndex, 	wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			//wxLogDebug(listControl->GetItemText(itemIndex));
			render_panel->raytracer.s.objects[obj_id]->remove_normal(firstSel);
		}
		//render_panel->raytracer.s.objects[obj_id]->remove_normal(item_id);
		render_panel->start_render();
		break;
	case ID_MOVEUP_NORMAL:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_normal(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_NORMAL:
		if (item_id >= render_panel->raytracer_app->m_NormalFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_normal(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDNULL_NORMAL:
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->add_null_normalmap();
		render_panel->start_render();
		break;
	case ID_CHANGE_COLOR_NORMAL:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_null_normalmap(item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_NORMAL:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_normalmap(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}
void RaytracerFrame::OnPopupClickAlpha(wxCommandEvent &evt) {
	int obj_id = render_panel->selected_object;
	if (obj_id < 0) return;
	if (obj_id >= render_panel->raytracer.s.objects.size()) return;
	int item_id = (int)static_cast<wxMenu *>(evt.GetEventObject())->GetClientData();
	int itemIndex = -1;
	int firstSel = render_panel->raytracer_app->m_AlphaFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
	switch (evt.GetId()) {
	case ID_ALBEDO_DELETE_ALPHA:
		render_panel->stop_render();
		while ((itemIndex = render_panel->raytracer_app->m_AlphaFile->GetNextItem(itemIndex, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) {
			render_panel->raytracer.s.objects[obj_id]->remove_alpha(firstSel);
		}
		render_panel->start_render();
		break;
	case ID_MOVEUP_ALPHA:
		if (item_id <= 0) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_alpha(item_id, item_id - 1);
		render_panel->start_render();
		break;
	case ID_MOVEDOWN_ALPHA:
		if (item_id >= render_panel->raytracer_app->m_AlphaFile->GetItemCount() - 1) return;
		render_panel->stop_render();
		render_panel->raytracer.s.objects[obj_id]->swap_alpha(item_id, item_id + 1);
		render_panel->start_render();
		break;
	case ID_ADDWHITE_ALPHA:
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->add_col_alpha(col.Red() / 255.);
			render_panel->start_render();		
		}
		break;
	case ID_CHANGE_COLOR_ALPHA:
	{
		if (render_panel->raytracer_app->colPicker->ShowModal() == wxID_OK) {
			wxColourData retData = render_panel->raytracer_app->colPicker->GetColourData();
			wxColour col = retData.GetColour();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_col_alpha(col.Red() / 255., item_id);
			render_panel->start_render();
		}
		break;
	}
	case ID_CHANGE_TEXTURE_ALPHA:
	{
		if (render_panel->raytracer_app->texOpenDlg->ShowModal() == wxID_OK) {
			std::string retData = render_panel->raytracer_app->texOpenDlg->GetPath().ToStdString();
			render_panel->stop_render();
			render_panel->raytracer.s.objects[obj_id]->set_alphamap(retData.c_str(), item_id);
			render_panel->start_render();
		}
		break;
	}
	}
	render_panel->update_gui();
}

void RaytracerFrame::OnListRightClick(wxListEvent &evt) {
	int sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP, "Move Up");
	mnu.Append(ID_MOVEDOWN, "Move Down");
	mnu.Append(ID_ALBEDO_DELETE, "Delete Row");
	mnu.Append(ID_REMOVE_TEXTURE, "Remove Texture");
	mnu.Append(ID_ADDWHITE, "Add Color");
	mnu.Append(ID_CHANGE_COLOR, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE, "Change Texture");
	
	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClick), NULL, this);
	PopupMenu(&mnu);
}

void RaytracerFrame::OnListRightClickSpecular(wxListEvent &evt) {
	int sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_SPECULAR, "Move Up");
	mnu.Append(ID_MOVEDOWN_SPECULAR, "Move Down");
	mnu.Append(ID_ALBEDO_DELETE_SPECULAR, "Delete Row");
	mnu.Append(ID_REMOVE_TEXTURE_SPECULAR, "Remove Texture");
	mnu.Append(ID_ADDWHITE_SPECULAR, "Add Color");
	mnu.Append(ID_CHANGE_COLOR_SPECULAR, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE_SPECULAR, "Change Texture");

	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickSpecular), NULL, this);
	PopupMenu(&mnu);
}
void RaytracerFrame::OnListRightClickNormal(wxListEvent &evt) {
	int sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_NORMAL, "Move Up");
	mnu.Append(ID_MOVEDOWN_NORMAL, "Move Down");
	mnu.Append(ID_ALBEDO_DELETE_NORMAL, "Delete");
	mnu.Append(ID_ADDNULL_NORMAL, "Add Null");
	mnu.Append(ID_CHANGE_COLOR_NORMAL, "Change to Null");
	mnu.Append(ID_CHANGE_TEXTURE_NORMAL, "Change Texture");

	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickNormal), NULL, this);
	PopupMenu(&mnu);
}
void RaytracerFrame::OnListRightClickAlpha(wxListEvent &evt) {
	int sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_ALPHA, "Move Up");
	mnu.Append(ID_MOVEDOWN_ALPHA, "Move Down");
	mnu.Append(ID_ALBEDO_DELETE_ALPHA, "Delete");
	mnu.Append(ID_ADDWHITE_ALPHA, "Add Color");
	mnu.Append(ID_CHANGE_COLOR_ALPHA, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE_ALPHA, "Change Texture");

	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickAlpha), NULL, this);
	PopupMenu(&mnu);
}
void RaytracerFrame::OnListRightClickRoughness(wxListEvent &evt) {
	int sel = evt.GetItem();
	void *data = reinterpret_cast<void *>(sel);
	wxMenu mnu;
	mnu.SetClientData(data);
	mnu.Append(ID_MOVEUP_ROUGHNESS, "Move Up");
	mnu.Append(ID_MOVEDOWN_ROUGHNESS, "Move Down");
	mnu.Append(ID_DELETE_ROUGHNESS, "Delete Row");
	mnu.Append(ID_REMOVE_TEXTURE_ROUGHNESS, "Remove Texture");
	mnu.Append(ID_ADDWHITE_ROUGHNESS, "Add Color");
	mnu.Append(ID_CHANGE_COLOR_ROUGHNESS, "Change Color");
	mnu.Append(ID_CHANGE_TEXTURE_ROUGHNESS, "Change Texture");
	mnu.Connect(wxEVT_COMMAND_MENU_SELECTED, wxCommandEventHandler(RaytracerFrame::OnPopupClickRoughness), NULL, this);
	PopupMenu(&mnu);
}

bool DnDFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_pOwner != NULL) {
		m_pOwner->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {		

			if (filenames[n].find(".xyz") != std::string::npos) {

				std::string values = wxGetTextFromUser("Space separated, -1: ignore, 0:x 1:y, 2:z, 3:nx,4:ny,5:nz, 6:colr, 7:colg, 8:colb",
					"Enter XYZ file format", "-1 -1 0 1 2 6 7 8");
				const char* line = &values.c_str()[0];
				int cols[100];
				int nbcols = 0;
				int offset;
				while (sscanf(line, "%d%n", &cols[nbcols], &offset) == 1) {
					nbcols++;
					line += offset;
				}

				PointSet* g = new PointSet(filenames[n], nbcols, cols, 0.005);
				g->scale = 30;
				Vector c = (g->bvh.bbox.bounds[0] + g->bvh.bbox.bounds[1])*0.5; // c is at 0
				double s = std::max(g->bvh.bbox.bounds[1][0] - g->bvh.bbox.bounds[0][0], std::max(g->bvh.bbox.bounds[1][1] - g->bvh.bbox.bounds[0][1], g->bvh.bbox.bounds[1][2] - g->bvh.bbox.bounds[0][2]));
				g->max_translation = Vector(0, m_pOwner->render_panel->raytracer.s.objects[2]->get_translation(1)[1] - (g->bvh.bbox.bounds[0][1])*g->scale, 0);
				m_pOwner->render_panel->raytracer.s.addObject(g);
			} else {
				Geometry* g = new Geometry(filenames[n], 1, Vector(0, 0, 0));
				g->scale = 30;
				g->display_edges = false;
				Vector c = (g->bvh.bbox.bounds[0] + g->bvh.bbox.bounds[1])*0.5; // c is at 0
				double s = std::max(g->bvh.bbox.bounds[1][0] - g->bvh.bbox.bounds[0][0], std::max(g->bvh.bbox.bounds[1][1] - g->bvh.bbox.bounds[0][1], g->bvh.bbox.bounds[1][2] - g->bvh.bbox.bounds[0][2]));
				g->max_translation = Vector(0, m_pOwner->render_panel->raytracer.s.objects[2]->get_translation(1)[1] - (g->bvh.bbox.bounds[0][1])*g->scale, 0);
				m_pOwner->render_panel->raytracer.s.addObject(g);
			}
		}		
		m_pOwner->render_panel->start_render();
	}
	return true;
}

bool DnDSpecularFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_specularmap(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}

bool DnDAlbedoFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_texture(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}

bool DnDNormalFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_normalmap(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}
bool DnDRoughnessFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_roughnessmap(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}
bool DnDAlphaFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (m_rtFrame->render_panel->selected_object < 0) return true;
	if (m_rtFrame->render_panel->selected_object >= m_rtFrame->render_panel->raytracer.s.objects.size()) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		for (size_t n = 0; n < nFiles; n++) {
			m_rtFrame->render_panel->raytracer.s.objects[m_rtFrame->render_panel->selected_object]->add_alphamap(filenames[n]);
		}
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}

bool DnDEnvmapFile::OnDropFiles(wxCoord, wxCoord, const wxArrayString& filenames)
{
	size_t nFiles = filenames.GetCount();
	if (nFiles < 1) return true;

	if (m_pOwner != NULL) {
		m_rtFrame->render_panel->stop_render();
		((Sphere*)m_rtFrame->render_panel->raytracer.s.objects[1])->load_envmap(filenames[0]);
		m_rtFrame->render_panel->update_gui();
		m_rtFrame->render_panel->start_render();
	}

	return true;
}

void RaytracerFrame::OnClose(wxCloseEvent& event)
{
	render_panel->stop_render();
	Destroy();  // you may also do:  event.Skip();
				// since the default event handler does call Destroy(), too
}