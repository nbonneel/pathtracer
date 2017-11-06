
#include "wx/wx.h"
#include "wx/dnd.h"
#include "wx/dirctrl.h"
#include "wx/notebook.h"
#include "wx/listctrl.h"
#include <wx/gauge.h>
#include <wx/clrpicker.h>
#include <wx/colordlg.h>
#include <wx/spinctrl.h>
#include <wx/button.h>
#include <ostream>
#include <sstream>
#include "chrono.h"
//#include <Windows.h>
#include <thread>

#include "Raytracer.h"

#define EDGES_CHECKBOX 1001
#define INTERP_CHECKBOX 1002
#define FOV_SLIDER 1003
#define APERTURE_SLIDER 1004
#define KS_SLIDER 1005
#define ALBEDO_FILES 1006
#define SPECULAR_FILES 10061
#define NORMAL_FILES 10062
#define ALPHA_FILES 10063
#define ROUGHNESS_FILES 10064
#define FILTER_SLIDER 1007
#define ALBEDO_COLORPICKER 1008
#define FOGDENSITY_SLIDER 1009
#define UNIFORMFOG_RADIO 1010
#define EXPFOG_RADIO 1010
#define ENVMAPINTENSITY_SLIDER 1011
#define LIGHTINTENSITY_SLIDER 1012
#define BOUNCES_SPIN 1013
#define TRANSPARENT_CHECKBOX 1014
#define FLIPNORMALS_CHECKBOX 1015
#define REFRACTION_INDEX 1016
#define DELETE_OBJECT 1017
#define RENDER_WIDTH 1018
#define RENDER_HEIGHT 1019
#define FOCUS_SLIDER 1020
#define NBRAYS 1021
#define LAUNCH_RENDER 1022

#define ID_ALBEDO_DELETE 10
#define ID_MOVEUP 11
#define ID_MOVEDOWN 12
#define ID_ADDWHITE 13
#define ID_ALBEDO_DELETE_ALPHA 14
#define ID_MOVEUP_ALPHA 15
#define ID_MOVEDOWN_ALPHA 16
#define ID_ADDWHITE_ALPHA 17
#define ID_ALBEDO_DELETE_NORMAL 18
#define ID_MOVEUP_NORMAL 19
#define ID_MOVEDOWN_NORMAL 20
#define ID_ADDNULL_NORMAL 21
#define ID_ALBEDO_DELETE_SPECULAR 22
#define ID_MOVEUP_SPECULAR 23
#define ID_MOVEDOWN_SPECULAR 24
#define ID_ADDWHITE_SPECULAR 25
#define ID_DELETE_ROUGHNESS 26
#define ID_MOVEUP_ROUGHNESS 27
#define ID_MOVEDOWN_ROUGHNESS 28
#define ID_ADDWHITE_ROUGHNESS 29
#define ID_SAVE 30
#define ID_OPEN 31
#define ID_CHANGE_COLOR 32
#define ID_CHANGE_TEXTURE 33
#define ID_CHANGE_COLOR_SPECULAR 34
#define ID_CHANGE_TEXTURE_SPECULAR 35
#define ID_CHANGE_COLOR_NORMAL 36
#define ID_CHANGE_TEXTURE_NORMAL 37
#define ID_CHANGE_COLOR_ALPHA 38
#define ID_CHANGE_TEXTURE_ALPHA 39
#define ID_CHANGE_COLOR_ROUGHNESS 40
#define ID_CHANGE_TEXTURE_ROUGHNESS 41
#define ID_REMOVE_TEXTURE 42
#define ID_REMOVE_TEXTURE_ROUGHNESS 43
#define ID_REMOVE_TEXTURE_SPECULAR 44
#define ID_SAVEIMAGE 45
#define ID_EXPORT_MTL 46


class RaytracerApp;
template<typename T> class RenderPanel;

class RaytracerFrame : public wxFrame
{
public:
	RaytracerFrame();
	virtual ~RaytracerFrame() {};


	void OnPaint(wxPaintEvent& event) {};
	void OnSize(wxSizeEvent& event) {
		Refresh();
		event.Skip();
	}
	void OnClose(wxCloseEvent& event);

	void OnLeftDown(wxMouseEvent& event) {};
	void OnRightDown(wxMouseEvent& event) {};

	void OnUpdateUIMoveByDefault(wxUpdateUIEvent& event) {};

	void OnPopupClick(wxCommandEvent &evt);
	void OnPopupClickSpecular(wxCommandEvent &evt);
	void OnPopupClickNormal(wxCommandEvent &evt);
	void OnPopupClickAlpha(wxCommandEvent &evt);
	void OnPopupClickRoughness(wxCommandEvent &evt);
	void OnListRightClick(wxListEvent &evt);
	void OnListRightClickSpecular(wxListEvent &evt);
	void OnListRightClickNormal(wxListEvent &evt);
	void OnListRightClickAlpha(wxListEvent &evt);
	void OnListRightClickRoughness(wxListEvent &evt);

	void SaveAs(wxCommandEvent &evt);
	void SaveImage(wxCommandEvent &evt);
	void ExportMtl(wxCommandEvent &evt);
	void Open(wxCommandEvent &evt);

	RenderPanel<double> *render_panel;
private:

	wxDECLARE_EVENT_TABLE();
};

class DnDFile : public wxFileDropTarget
{
public:
    DnDFile(RaytracerFrame *pOwner = NULL) { m_pOwner = pOwner; }

    virtual bool OnDropFiles(wxCoord x, wxCoord y,
                             const wxArrayString& filenames) wxOVERRIDE;

private:
	RaytracerFrame *m_pOwner;
};

class DnDAlbedoFile : public wxFileDropTarget
{
public:
	DnDAlbedoFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) wxOVERRIDE;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDSpecularFile : public wxFileDropTarget
{
public:
	DnDSpecularFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) wxOVERRIDE;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDRoughnessFile : public wxFileDropTarget
{
public:
	DnDRoughnessFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) wxOVERRIDE;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDEnvmapFile : public wxFileDropTarget
{
public:
	DnDEnvmapFile(wxTextCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) wxOVERRIDE;

private:
	wxTextCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDNormalFile : public wxFileDropTarget
{
public:
	DnDNormalFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) wxOVERRIDE;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

class DnDAlphaFile : public wxFileDropTarget
{
public:
	DnDAlphaFile(wxListCtrl* pOwner = NULL, RaytracerFrame *rtFrame = NULL) {
		m_pOwner = pOwner; m_rtFrame = rtFrame;
	}

	virtual bool OnDropFiles(wxCoord x, wxCoord y,
		const wxArrayString& filenames) wxOVERRIDE;

private:
	wxListCtrl *m_pOwner;
	RaytracerFrame *m_rtFrame;
};

wxBEGIN_EVENT_TABLE(RaytracerFrame, wxFrame)
EVT_LEFT_DOWN(RaytracerFrame::OnLeftDown)
EVT_RIGHT_DOWN(RaytracerFrame::OnRightDown)
EVT_PAINT(RaytracerFrame::OnPaint)
EVT_CLOSE(RaytracerFrame::OnClose)
EVT_SIZE(RaytracerFrame::OnSize)
wxEND_EVENT_TABLE()



static void doRender(void* Param);



template<typename T>
class RenderPanel : public wxPanel
{

public:

	RenderPanel(RaytracerFrame* parent, RaytracerApp* app):wxPanel(parent) {
		Connect(wxEVT_PAINT, wxPaintEventHandler(RenderPanel::paintEvent));

		raytracer_app = app;
		should_terminate = false;
		raytracer.loadScene();
		W = raytracer.W;
		H = raytracer.H;
		cur_img.resize(W*H * 3);
		for (int i = 0; i < H; i++) {
			for (int j = 0; j < W; j++) {
				cur_img[(i*W + j) * 3 + 0] = 255;
				cur_img[(i*W + j) * 3 + 1] = 0;
				cur_img[(i*W + j) * 3 + 2] = 0;
			}
		}

		left_mouse_down = false;
		parent->render_panel = this;
		wasDragging = false;

		//compute_thread = CreateThread(NULL, 0, doRender, (void*) this, 0, NULL);
		//SetThreadPriority(compute_thread, THREAD_PRIORITY_HIGHEST);
		start_render();
		//SetThreadPriority(compute_thread, THREAD_PRIORITY_LOWEST);
	}


	/*void partition(int* indices, int i0, int i1) {

		int mini = 100, maxi = -100;
		for (int i = i0; i < i1; i++) {
			mini = std::min(mini, indices[i]);
			maxi = std::max(maxi, indices[i]);
		}
		double split_val = (mini + maxi) / 2;

		int pivot = i0 - 1;

		for (int i = i0; i < i1; i++) {
			double center_split_dim = indices[i];

			if (center_split_dim <= split_val) {
				pivot++;
				std::swap(indices[i], indices[pivot]);
			}
		}

		if (pivot < i0 || pivot >= i1 - 1 || i1 <= i0 + 1) {
			return;
		}

		partition(indices, i0, pivot + 1);
		partition(indices, pivot + 1, i1);

	}


	int test() {

		int values[12];
		for (int k = 0; k < 20; k++) {
			for (int i = 0; i < 12; i++) {
				values[i] = rand() % 100;
			}
			partition(values, 0, 12);
		}
		return 0;
	}*/

	void paintEvent(wxPaintEvent& evt)
	{
	   //test();
		wxPaintDC dc(this);
		render(dc);
	}

	inline double fastPow(double a, double b) {
		union {
			double d;
			int x[2];
		} u = { a };
		u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
		u.x[0] = 0;
		return u.d;
	}

	void delete_object(wxCommandEvent& event) {
		if (selected_object < 0) return;
		if (selected_object >= raytracer.s.objects.size()) return;

		stop_render();
		raytracer.s.deleteObject(selected_object);
		start_render();
	}

	void launch_render(wxCommandEvent& event) {
		stop_render();
		PerfChrono perf;
		perf.Start();
		raytracer.stopped = false;
		raytracer.render_image();
		raytracer.stopped = true;
		double t = perf.GetDiffMs();
		FILE* f = fopen("time.txt", "a+");
		fprintf(f, "%e\n", t);
		fclose(f);
		start_render();
	}

	void update_parameters_and_render(wxCommandEvent& event) {
		
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
		raytracer.cam.aperture = raytracer_app->aperture_slider->GetValue()/1000.;
		raytracer.cam.focus_distance = raytracer_app->focus_slider->GetValue() / 100.;
		raytracer.sigma_filter = raytracer_app->filter_slider->GetValue() / 10.;
		//wxColour alb = raytracer_app->albedoColorPicker->GetColour();
		//raytracer.s.objects[selected_object]->albedo = Vector(alb.Red()/255., alb.Green()/255., alb.Blue()/255.);

		//raytracer.s.objects[selected_object]->ks = raytracer_app->ks_slider->GetValue() / 100.;
		raytracer.s.fog_density = raytracer_app->fogdensity_slider->GetValue() / 100.;

		raytracer.s.intensite_lumiere = raytracer_app->lightintensity_slider->GetValue() / 100. * 1000000000 * 4.*M_PI / (4.*M_PI*raytracer.s.lumiere->R*raytracer.s.lumiere->R*M_PI);
		raytracer.s.envmap_intensity = raytracer_app->envmapintensity_slider->GetValue()/100.;

		raytracer.s.fog_type = raytracer_app->uniformFogRadio->GetValue()?0:1;
		raytracer.nb_bounces = raytracer_app->bounces->GetValue();
		
		raytracer.s.objects[selected_object]->flip_normals = raytracer_app->flipnormals->IsChecked();
		raytracer.s.objects[selected_object]->refr_index = raytracer_app->refractionIndex->GetValue();

		

		start_render();
	}
	void update_textures_and_render(wxCommandEvent& event) {

		if (selected_object < 0) return;
		if (selected_object >= raytracer.s.objects.size()) return;

		stop_render();

		raytracer.s.objects[selected_object]->display_edges = raytracer_app->show_edges->IsChecked();
		
		start_render();
	}

	void update_gui() {

		if ((selected_object < 0) || (selected_object >= raytracer.s.objects.size())) {
			raytracer_app->objectName->SetLabelText(" ");
			return;
		}

		Geometry* g = dynamic_cast<Geometry*>(raytracer.s.objects[selected_object]);
		if (g) {
			std::ostringstream os;
			os << "triangles: " << g->indices.size() << ", vertices: " << g->vertices.size() <<", bvh leaves size: "<<g->max_bvh_triangles<<", bvh max depth: "<<g->bvh_depth<<", bvh avg depth:"<<g->bvh_avg_depth<<", bvh nb nodes: "<<g->bvh_nb_nodes<<", mouse distance: "<< selected_object_t;
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
		raytracer_app->show_edges->SetValue( raytracer.s.objects[selected_object]->display_edges );
		raytracer_app->interp_normals->SetValue(raytracer.s.objects[selected_object]->interp_normals);
		raytracer_app->fov_slider->SetValue(raytracer.cam.fov * 180. / M_PI);
		raytracer_app->aperture_slider->SetValue(raytracer.cam.aperture*1000);
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
		raytracer_app->envmapName->SetLabelText( ((Sphere*)raytracer.s.objects[1])->envmapfilename);

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
		raytracer_app->uniformFogRadio->SetValue(raytracer.s.fog_type==0);

		double factor = 1. / 100. * 1000000000 * 4.*M_PI / (4.*M_PI*raytracer.s.lumiere->R*raytracer.s.lumiere->R*M_PI);
		raytracer_app->lightintensity_slider->SetValue(raytracer.s.intensite_lumiere / factor);
		raytracer_app->envmapintensity_slider->SetValue(raytracer.s.envmap_intensity * 100);
	}



	std::vector<double> extrapolated_image;
	std::vector<bool> computed, computed2;
	void render(wxDC& dc)
	{
		if (cur_img.size() == 0) return;
		if (raytracer.stopped) return;
		raytracer_app->progressBar->SetValue(raytracer.current_nb_rays/ (double)raytracer.nrays*1000.);
		std::ostringstream os;
		os << "Time per ray: " << raytracer.curTimePerFrame/1000. <<" s";
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
				if (raytracer.sample_count[i*raytracer.W + j]<=5) {
					double alpha = raytracer.sample_count[i*raytracer.W + j] / 6.;

					int ix = j / 16;
					double fx = (j / 16.) - ix;

					int iy = i / 16;
					double fy = (i / 16.) - iy;

					for (int k = 0; k < 3; k++) {
						//double lowval = raytracer.imagedouble_lowres[(iy * raytracer.Wlr + ix) * 3 + k];  // nearest neighbor interp.
						double lowval;
						if (ix < raytracer.Wlr-1 && iy < raytracer.Hlr-1)
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
				double normalization = 1./(raytracer.sample_count[i*raytracer.W + j]+1.);
				/*gamma_corrected_image[(i*raytracer.W + j) * 3 + 0] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 0] / (raytracer.current_nb_rays + 1)), 1 / 2.2));   // rouge
				gamma_corrected_image[(i*raytracer.W + j) * 3 + 1] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 1] / (raytracer.current_nb_rays + 1)), 1 / 2.2)); // vert
				gamma_corrected_image[(i*raytracer.W + j) * 3 + 2] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 2] / (raytracer.current_nb_rays + 1)), 1 / 2.2)); // bleu*/
				gamma_corrected_image[(i*raytracer.W + j) * 3 + 0] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 0] * normalization), 1 / 2.2));   // rouge
				gamma_corrected_image[(i*raytracer.W + j) * 3 + 1] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 1] * normalization), 1 / 2.2)); // vert
				gamma_corrected_image[(i*raytracer.W + j) * 3 + 2] = std::min(255., fastPow(std::max(0., extrapolated_image[(i*raytracer.W + j) * 3 + 2] * normalization), 1 / 2.2)); // bleu

			}
		}

		static int nbcalls = -1; nbcalls++;
		if (nbcalls == 0 || screenImage.GetWidth()!= raytracer.W || screenImage.GetHeight()!= raytracer.H) {
			screenImage = wxImage(raytracer.W, raytracer.H, &(gamma_corrected_image[0]), true);
		}

		bmpBuf = wxBitmap(screenImage, dc);
		

		double scale_x = (double)displayW / raytracer.W;
		double scale_y = (double)displayH / raytracer.H;
		dc.SetUserScale(scale_x, scale_y);
		dc.DrawBitmap(bmpBuf, 0,0);
		dc.SetUserScale(1.0, 1.0);
	}



	void paintNow()
	{
		wxClientDC dc(this);
		render(dc);
	}
	void stop_render() {
		raytracer.stopRender();		
		//WaitForSingleObject(compute_thread, INFINITE);
		compute_thread2.join();
	}

	void start_render() {
		//compute_thread = CreateThread(NULL, 0, doRender, (void*)this, 0, NULL);
		//SetThreadPriority(compute_thread, THREAD_PRIORITY_HIGHEST);
		raytracer.clear_image();
		raytracer.stopped = false;
		compute_thread2 = std::thread(doRender, (void*)this);		
	}

	void mouseDown(wxMouseEvent& event) {
		mouse_init_x = event.m_x; // coordinates before dragging the mouse
		mouse_init_y = event.m_y;
	}

	void mouseMoved(wxMouseEvent& event)
	{
		mouse_x = event.m_x;
		mouse_y = event.m_y;

		if (mouse_x > displayW + 10) return;
		if (mouse_y > displayH + 10) return;

		if ((event.LeftIsDown() || event.RightIsDown() || event.MiddleIsDown()) && event.Dragging())
		{
			left_mouse_down = event.LeftIsDown();
			right_mouse_down = event.RightIsDown();
			middle_mouse_down = event.MiddleIsDown();
			alt_down = event.AltDown();
			shift_down = event.ShiftDown();
			wasDragging = true;

			int dx = mouse_x - mouse_prev_x;
			int dy = -mouse_y + mouse_prev_y;
			int Dx = mouse_x - mouse_init_x;
			int Dy = -mouse_y + mouse_init_y;
			
			if (shift_down) {
				if (left_mouse_down) {
					if (selected_object >= 0 && selected_object < raytracer.s.objects.size()) {
						//raytracer.s.objects[selected_object]->max_translation += Vector(dx, dy, 0.)*0.3;
						Vector camera_right = cross(raytracer.cam.direction, raytracer.cam.up);
						raytracer.s.objects[selected_object]->max_translation += dx*selected_object_t*tan(raytracer.cam.fov/ displayW)*camera_right + dy*selected_object_t*tan(raytracer.cam.fov / displayH)*raytracer.cam.up;
					}
				}
				if (middle_mouse_down) {
					if (selected_object >= 0 && selected_object < raytracer.s.objects.size()) {
						raytracer.s.objects[selected_object]->max_rotation += Vector(0., 0., dy*M_PI / 180.*3.);
					}
				}
				if (right_mouse_down) {
					if (selected_object >= 0 && selected_object < raytracer.s.objects.size()) {
						raytracer.s.objects[selected_object]->max_rotation += Vector(dy*M_PI / 180.*3., dx*M_PI / 180.*3., 0.);
					}
				}
			} else { // camera motion
				if (left_mouse_down) {

					if (abs(mouse_x - displayW) < 10 && abs(mouse_y - displayH) < 10) {
						displayW = mouse_x;
						displayH = mouse_y;
					} else
						raytracer.cam.rotateAroundAxes(dx*M_PI / 180 * 1, dy*M_PI / 180 * 1, 1);
				}
				if (right_mouse_down) {
					Vector camera_right = cross(raytracer.cam.direction, raytracer.cam.up);
					raytracer.cam.position += dx*camera_right + dy*raytracer.cam.up;

				}
			}			
			stop_render();
			start_render();
										
		} else {
			if (!event.Dragging()) {

			}
		}
		mouse_prev_x = mouse_x;
		mouse_prev_y = mouse_y;

		paintNow();
	}

	void mouseUp(wxMouseEvent& event) {
				
		if ( (event.m_x > displayW + 10) || (event.m_y > displayW + 10)) {
			wasDragging = false;
			return;
		}	
		

		/*if (!wasDragging)*/ {
			Ray r = raytracer.cam.generateDirection((displayH - (mouse_init_y - 1.))*(double)raytracer.H / displayH, mouse_init_x*(double)raytracer.W / displayW, 0, 0, 0, 0, 0, raytracer.W, raytracer.H);

			Vector P;
			int new_selected, new_tri;
			double new_t;
			MaterialValues newmat;
			bool has_inter = raytracer.s.intersection(r, P, new_selected, new_t, newmat, new_tri);
			if (has_inter) {
				selected_object = new_selected;
				selected_object_t = new_t;
				selected_tri = new_tri;
				update_gui();
			}			
			mouse_prev_x = mouse_x;
			mouse_prev_y = mouse_y;

		}
		
		left_mouse_down = false;
		if (!event.AltDown())
			alt_down = false;

		wasDragging = false;
	}

	void mouseWheel(wxMouseEvent& event) {
		if (event.ShiftDown()) {
			if (selected_object >= 0 && selected_object < raytracer.s.objects.size()) {
				raytracer.s.objects[selected_object]->scale *= (event.GetWheelRotation() > 0) ? 1.1 : (1 / 1.1);
				stop_render();
				start_render();
			}
		} else {
			double dist =  (event.GetWheelRotation() > 0 ? 1 : -1)*2.;
			raytracer.cam.position += dist*raytracer.cam.direction;
			stop_render();
			start_render();
		}
		paintNow();
	}

	bool should_terminate;
	bool left_mouse_down, right_mouse_down, middle_mouse_down, alt_down, shift_down, wasDragging;
	int mouse_x, mouse_y, mouse_init_x, mouse_init_y, mouse_prev_x, mouse_prev_y, selected_object, selected_tri;
	double selected_object_t;
	std::vector<unsigned char> cur_img;
	int W, H;
	int displayW, displayH;
	//HANDLE compute_thread;
	std::thread compute_thread2;

	RaytracerApp* raytracer_app;
	Raytracer raytracer;
	wxImage screenImage;
	wxBitmap bmpBuf;
	std::vector<unsigned char> gamma_corrected_image;

	DECLARE_EVENT_TABLE()
};

static void doRender(void* Param) {
	Raytracer* raytracer = &(((RenderPanel<double>*)Param)->raytracer);
	raytracer->render_image();
}

BEGIN_EVENT_TABLE_TEMPLATE1(RenderPanel, wxPanel, T)
EVT_PAINT(RenderPanel<T>::paintEvent)
EVT_MOTION(RenderPanel<T>::mouseMoved)
EVT_LEFT_DOWN(RenderPanel<T>::mouseDown)
EVT_RIGHT_DOWN(RenderPanel<T>::mouseDown)
EVT_LEFT_UP(RenderPanel<T>::mouseUp)
EVT_RIGHT_UP(RenderPanel<T>::mouseUp)
EVT_MOUSEWHEEL(RenderPanel<T>::mouseWheel)
END_EVENT_TABLE()


class RaytracerApp : public wxApp
{
public:
	virtual bool OnInit() wxOVERRIDE;
	RenderPanel<double>* renderPanel;
	wxGauge *progressBar;
	wxBookCtrlBase *m_bookCtrl;
	wxCheckBox *show_edges, *interp_normals, *transparent, *flipnormals;
	wxTextCtrl *objectName, *envmapName;
	wxSlider *fov_slider, *aperture_slider, /**ks_slider,*/ *filter_slider, *fogdensity_slider, *envmapintensity_slider, *lightintensity_slider, *focus_slider;
	wxListCtrl *m_AlbedoFile, *m_SpecularFile, *m_NormalFile, *m_AlphaFile, *m_RoughnessFile;
	//wxColourPickerCtrl *albedoColorPicker;
	wxRadioButton *uniformFogRadio, *expFogRadio;
	wxColourDialog* colPicker;
	wxFileDialog* texOpenDlg;
	wxSpinCtrl *bounces, *renderwidth, *renderheight, *nbrays;
	wxSpinCtrlDouble* refractionIndex;
	wxButton *deleteObject, *launchRender;
	wxStaticText* infoModel;
	wxStaticText* infoPerf;
	bool render_loop_on;

	void activateRenderLoop(bool on)
	{
		renderPanel->SetDoubleBuffered(true);
		if (on && !render_loop_on)
		{
			Connect(wxID_ANY, wxEVT_IDLE, wxIdleEventHandler(RaytracerApp::OnIdle));
			Connect(wxID_ANY, wxEVT_ERASE_BACKGROUND, wxEraseEventHandler(RaytracerApp::OnEraseBackGround));
			render_loop_on = true;
		} else if (!on && render_loop_on)
		{
			Disconnect(wxEVT_IDLE, wxIdleEventHandler(RaytracerApp::OnIdle));
			render_loop_on = false;
		}
	}

	void OnEraseBackGround(wxEraseEvent& event) {
		//wxDC * TheDC = event.GetDC();		
		event.Skip(); // don't call -- not needed, we already erased the background
	}

	void OnIdle(wxIdleEvent& evt)
	{
		if (render_loop_on)
		{
			renderPanel->paintNow();
			evt.RequestMore(); // render continuously, not only once on idle				
		}
	}
	void OnUnhandledException() {
		throw;
	};

	DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(RaytracerApp, wxApp)
EVT_IDLE(OnIdle)
END_EVENT_TABLE()



