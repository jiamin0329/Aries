


/*!
* \struct GEOM_GeometryMultigridQueue
* \brief Class for a multigrid queue system
* \author F. Palacios
* \version 3.2.9 "eagle"
* \date Aug 12, 2012
*/


#ifndef ARIES_GEOM_GEOMETRYMULTIGRIDQUEUE_HPP
#define ARIES_GEOM_GEOMETRYMULTIGRIDQUEUE_HPP

#include "../Common/TBOX_Config.hpp"
#include "GEOM_Geometry.hpp"

namespace ARIES
{
    namespace GEOM
    {
        class GEOM_GeometryMultigridQueue 
        {
            std::vector<std::vector<unsigned long> > QueueCV; /*!< \brief Queue structure to choose the next control volume in the agglomeration process. */
            short *Priority;	/*!< \brief The priority is based on the number of pre-agglomerated neighbors. */
            bool *RightCV;	/*!< \brief In the lowest priority there are some CV that can not be agglomerated, this is the way to identify them */
            unsigned long nPoint; /*!< \brief Total number of points. */

        public:

            /*!
            * \brief Constructor of the class.
            * \param[in] val_npoint - Number of control volumes.
            */
            GEOM_GeometryMultigridQueue(unsigned long val_npoint);

            /*!
            * \brief Destructor of the class.
            */
            ~GEOM_GeometryMultigridQueue(void);

            /*!
            * \brief Add a new CV to the list.
            * \param[in] val_new_point - Index of the new point.
            * \param[in] val_number_neighbors - Number of neighbors of the new point.
            */
            void AddCV(unsigned long val_new_point, unsigned short val_number_neighbors);

            /*!
            * \brief Remove a CV from the list.
            * \param[in] val_remove_point - Index of the control volume to be removed.
            */
            void RemoveCV(unsigned long val_remove_point);

            /*!
            * \brief Change a CV from a list to a different list.
            * \param[in] val_move_point - Index of the control volume to be moved.
            * \param[in] val_number_neighbors - New number of neighbors of the control volume.
            */
            void MoveCV(unsigned long val_move_point, short val_number_neighbors);

            /*!
            * \brief Increase the priority of the CV.
            * \param[in] val_incr_point - Index of the control volume.
            */
            void IncrPriorityCV(unsigned long val_incr_point);

            /*!
            * \brief Increase the priority of the CV.
            * \param[in] val_red_point - Index of the control volume.
            */
            void RedPriorityCV(unsigned long val_red_point);

            /*!
            * \brief Visualize the control volume queue.
            */
            void VisualizeQueue(void);

            /*!
            * \brief Visualize the priority list.
            */
            void VisualizePriority(void);

            /*!
            * \brief Find a new seed control volume.
            * \return Index of the new control volume.
            */
            long NextCV(void);

            /*!
            * \brief Check if the queue is empty.
            * \return <code>TRUE</code> or <code>FALSE</code> depending if the queue is empty.
            */
            bool EmptyQueue(void);

            /*!
            * \brief Total number of control volume in the queue.
            * \return Total number of control points.
            */
            unsigned long TotalCV(void);

            /*!
            * \brief Update the queue with the new control volume (remove the CV and
            increase the priority of the neighbors).
            * \param[in] val_update_point - Index of the new point.
            * \param[in] fine_grid - Fine grid geometry.
            */
            void Update(unsigned long val_update_point, GEOM_Geometry *fine_grid);
        };
    }
}

#include "GEOM_GeometryMultigridQueue.inl"

#endif