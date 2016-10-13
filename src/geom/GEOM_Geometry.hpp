



namespace ARIES
{
    namespace GEOM
    {
        class GEOM_Geometry
        {
        public:
            /*!
             * \brief Constructor of the class.
             */
            GEOM_Geometry(void);

            /*!
             * \brief Destructor of the class.
             */
            virtual ~GEOM_Geometry(void);

            unsigned long GetnLine(void);




            /*!
             * \brief Get number of vertices.
             * \param[in] val_marker - Marker of the boundary.
             * \return Number of vertices.
             */
            unsigned long GetnVertex(unsigned short val_marker);

            /*!
             * \brief Get the edge index from using the nodes of the edge.
             * \param[in] first_point - First point of the edge.
             * \param[in] second_point - Second point of the edge.
             * \return Index of the edge.
             */
            long FindEdge(unsigned long first_point, unsigned long second_point);

            /*!
             * \brief Get the edge index from using the nodes of the edge.
             * \param[in] first_point - First point of the edge.
             * \param[in] second_point - Second point of the edge.
             * \return Index of the edge.
             */
            bool CheckEdge(unsigned long first_point, unsigned long second_point);

            
            /*!
             * \brief Create a file for testing the geometry.
             */
            void TestGeometry(void);


            /*!
             * \brief Get the index of a marker.
             * \param[in] val_marker - Marker of the boundary.
             * \return Index of the marker in the grid defintion.
             */
            std::string GetMarker_Tag(unsigned short val_marker);

            /*!
             * \brief Set index of a marker.
             * \param[in] val_marker - Marker of the boundary.
             * \param[in] val_index - Index of the marker.
             */
            void SetMarker_Tag(unsigned short val_marker, std::string val_index);


  
            /*!
             * \brief A virtual function.
             * \param[in] first_elem - Identification of the first element.
             * \param[in] second_elem - Identification of the second element.
             * \param[in] face_first_elem - Index of the common face for the first element.
             * \param[in] face_second_elem - Index of the common face for the second element.
             */
            virtual bool FindFace(unsigned long first_elem, unsigned long second_elem, unsigned short &face_first_elem,
                unsigned short &face_second_elem);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void ComputeWall_Distance(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetPositive_ZArea(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             */
            virtual void SetPoint_Connectivity(void);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetRCM_Ordering(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             */
            virtual void SetElement_Connectivity(void);


            /*!
             * \brief A virtual member.
             */
            virtual void SetBoundVolume(void);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetVertex(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             */
            virtual void SetVertex(void);

            /*!
             * \brief A virtual member.
             */
            virtual void SetCG(void);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] action - Allocate or not the new elements.
             */
            virtual void SetControlVolume(TBOX::TBOX_Config *config, unsigned short action);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] action - Allocate or not the new elements.
             */
            virtual void VisualizeControlVolume(TBOX::TBOX_Config *config, unsigned short action);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void MatchNearField(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void MatchActuator_Disk(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void MatchInterface(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] geometry_donor - Geometry of the donor zone.
             * \param[in] config_donor - Definition of the donor problem.
             */
            virtual void MatchZone(TBOX::TBOX_Config *config, GEOM_Geometry *geometry_donor, TBOX::TBOX_Config *config_donor,
                unsigned short val_iZone, unsigned short val_nZone);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] action - Allocate or not the new elements.
             */
            virtual void SetBoundControlVolume(TBOX::TBOX_Config *config, unsigned short action);

            /*!
             * \brief A virtual member.
             * \param[in] config_filename - Name of the file where the tecplot information is going to be stored.
             */
            virtual void SetTecPlot(char config_filename[TBOX::MAX_STRING_SIZE], bool new_file);

            /*!
             * \brief A virtual member.
             * \param[in] mesh_filename - Name of the file where the tecplot information is going to be stored.
             * \param[in] new_file - Boolean to decide if aopen a new file or add to a old one
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetBoundTecPlot(char mesh_filename[TBOX::MAX_STRING_SIZE], bool new_file, TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] mesh_filename - Name of the file where the tecplot information is going to be stored.
             * \param[in] new_file - Boolean to decide if aopen a new file or add to a old one
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetBoundSTL(char mesh_filename[TBOX::MAX_STRING_SIZE], bool new_file, TBOX::TBOX_Config *config);


            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void Check_IntElem_Orientation(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void Check_BoundElem_Orientation(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetColorGrid(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetColorGrid_Parallel(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void DivideConnectivity(TBOX::TBOX_Config *config, unsigned short Elem_Type);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetPeriodicBoundary(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_domain - Number of domains for parallelization purposes.
             */
            virtual void SetSendReceive(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_domain - Number of domains for parallelization purposes.
             */
            virtual void SetBoundaries(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] geometry - Geometrical definition of the problem.
             */
            virtual void SetCoord(GEOM_Geometry *geometry);

            /*!
             * \brief A virtual member.
             * \param[in] val_nSmooth - Number of smoothing iterations.
             * \param[in] val_smooth_coeff - Relaxation factor.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetCoord_Smoothing(unsigned short val_nSmooth, double val_smooth_coeff, TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] geometry - Geometrical definition of the problem.
             */
            virtual void SetPoint_Connectivity(GEOM_Geometry *geometry);

            /*!
             * \brief A virtual member.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetVertex(GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] action - Allocate or not the new elements.
             */
            virtual void SetControlVolume(TBOX::TBOX_Config *config, GEOM_Geometry *geometry, unsigned short action);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] action - Allocate or not the new elements.
             */
            virtual void SetBoundControlVolume(TBOX::TBOX_Config *config, GEOM_Geometry *geometry, unsigned short action);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_mesh_out_filename - Name of the output file.
             */
            virtual void SetMeshFile(TBOX::TBOX_Config *config, std::string val_mesh_out_filename);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] val_mesh_out_filename - Name of the output file.
             */
            virtual void SetMeshFile(GEOM_Geometry *geometry, TBOX::TBOX_Config *config, std::string val_mesh_out_filename);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetBoundSensitivity(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] geometry - Geometrical definition of the problem.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetPeriodicBoundary(GEOM_Geometry *geometry, TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetRotationalVelocity(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] iter - Current physical time step.
             */
            virtual void SetGridVelocity(TBOX::TBOX_Config *config, unsigned long iter);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void Set_MPI_Coord(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void Set_MPI_GridVel(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] geometry - Geometry of the fine mesh.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void SetRestricted_GridVelocity(GEOM_Geometry *fine_mesh, TBOX::TBOX_Config *config);

            /*!
             * \brief Find and store all vertices on a sharp corner in the geometry.
             * \param[in] config - Definition of the particular problem.
             */
            void ComputeSurf_Curvature(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            void ComputeAirfoil_Section(double *Plane_P0, double *Plane_Normal,
                double MinXCoord, double MaxXCoord, double *FlowVariable,
                std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil,
                std::vector<double> &Zcoord_Airfoil, std::vector<double> &Variable_Airfoil,
                bool original_surface, TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual double Compute_MaxThickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, TBOX::TBOX_Config *config, std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil, std::vector<double> &Zcoord_Airfoil, bool original_surface);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual double Compute_AoA(double *Plane_P0, double *Plane_Normal, unsigned short iSection, std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil, std::vector<double> &Zcoord_Airfoil, bool original_surface);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual double Compute_Chord(double *Plane_P0, double *Plane_Normal, unsigned short iSection, std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil, std::vector<double> &Zcoord_Airfoil, bool original_surface);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
             * \returns The minimum value of the airfoil thickness.
             */
            virtual double Compute_Thickness(double *Plane_P0, double *Plane_Normal, unsigned short iSection, double Location, TBOX::TBOX_Config *config, std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil, std::vector<double> &Zcoord_Airfoil, bool original_surface);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
             * \returns The total volume of the airfoil.
             */
            virtual double Compute_Area(double *Plane_P0, double *Plane_Normal, unsigned short iSection, TBOX::TBOX_Config *config, std::vector<double> &Xcoord_Airfoil, std::vector<double> &Ycoord_Airfoil, std::vector<double> &Zcoord_Airfoil, bool original_surface);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             * \param[in] original_surface - <code>TRUE</code> if this is the undeformed surface; otherwise <code>FALSE</code>.
             * \returns The total volume of the 3D body.
             */
            virtual double Compute_Volume(TBOX::TBOX_Config *config, bool original_surface);

            /*!
             * \brief A virtual member.
             * \param[in] config - Definition of the particular problem.
             */
            virtual void FindNormal_Neighbor(TBOX::TBOX_Config *config);

            /*!
             * \brief A virtual member.
             * \param[in] val_ipoint - Global point.
             * \returns Local index that correspond with the global index.
             */
            virtual long GetGlobal_to_Local_Point(long val_ipoint);

            /*!
             * \brief A virtual member.
             * \param[in] val_ipoint - Global marker.
             * \returns Local marker that correspond with the global index.
             */
            virtual unsigned short GetGlobal_to_Local_Marker(unsigned short val_imarker);




