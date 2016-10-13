




#include "GEOM_GeometryMultigridQueue.hpp"

namespace ARIES
{
    namespace GEOM
    {
        GEOM_GeometryMultigridQueue::GEOM_GeometryMultigridQueue(unsigned long val_npoint) 
        {
            unsigned long iPoint;

            nPoint = val_npoint;
            Priority = new short[nPoint];
            RightCV = new bool[nPoint];

            QueueCV.resize(1);

            /*--- Queue initialization with all the points in the finer grid ---*/
            for (iPoint = 0; iPoint < nPoint; iPoint++) 
            {
                QueueCV[0].push_back(iPoint);
                Priority[iPoint] = 0;
                RightCV[iPoint] = true;
            }
        }

        GEOM_GeometryMultigridQueue::~GEOM_GeometryMultigridQueue(void) 
        {
            delete[] Priority;
            delete[] RightCV;
        }

        void GEOM_GeometryMultigridQueue::AddCV(unsigned long val_new_point, unsigned short val_number_neighbors) 
        {
            unsigned short Max_Neighbors = QueueCV.size() - 1;

            /*--- Basic check ---*/
            if (val_new_point > nPoint)
            {
                std::cout << "The index of the CV is greater than the size of the priority list." << std::endl;
                exit(EXIT_FAILURE);
            }

            /*--- Resize the list ---*/
            if (val_number_neighbors > Max_Neighbors)
                QueueCV.resize(val_number_neighbors + 1);

            /*--- Find the point in the queue ---*/
            bool InQueue = false;
            if (Priority[val_new_point] == val_number_neighbors) InQueue = true;

            if (!InQueue) 
            {
                /*--- Add the control volume, and update the priority list ---*/
                QueueCV[val_number_neighbors].push_back(val_new_point);
                Priority[val_new_point] = val_number_neighbors;
            }
        }

        void GEOM_GeometryMultigridQueue::RemoveCV(unsigned long val_remove_point) 
        {
            unsigned short iPoint;
            bool check;

            /*--- Basic check ---*/
            if (val_remove_point > nPoint) 
            {
                std::cout << "The index of the CV is greater than the size of the priority list." << std::endl;
                exit(EXIT_FAILURE);
            }

            /*--- Find priority of the Control Volume ---*/
            short Number_Neighbors = Priority[val_remove_point];
            if (Number_Neighbors == -1) 
            {
                std::cout << "The CV " << val_remove_point << " is not in the priority list. (RemoveCV)" << std::endl;
                exit(EXIT_FAILURE);
            }

            /*--- Find the point in the queue ---*/
            std::vector<unsigned long>::iterator ItQueue = find(QueueCV[Number_Neighbors].begin(),
                QueueCV[Number_Neighbors].end(),
                val_remove_point);
            if (ItQueue != QueueCV[Number_Neighbors].end()) QueueCV[Number_Neighbors].erase(ItQueue);

            Priority[val_remove_point] = -1;

            /*--- Check that the size of the queue is the right one ---*/
            unsigned short Size_QueueCV = 0;
            check = false;
            for (iPoint = 0; iPoint < QueueCV.size(); iPoint++)
                if (QueueCV[iPoint].size() != 0) 
                { 
                    Size_QueueCV = iPoint; check = true; 
                }

            /*--- Resize the queue, if check = false, the queue is empty, at least
            we need one element in the queue ---*/
            if (check) QueueCV.resize(Size_QueueCV + 1);
            else QueueCV.resize(1);
        }

        void GEOM_GeometryMultigridQueue::MoveCV(unsigned long val_move_point, short val_number_neighbors) 
        {
            if (val_number_neighbors < 0) 
            {
                val_number_neighbors = 0;
                RightCV[val_move_point] = false;
            }
            else 
            {
                RightCV[val_move_point] = true;
            }

            /*--- Remove the control volume ---*/
            RemoveCV(val_move_point);

            /*--- Add a new control volume ---*/
            AddCV(val_move_point, val_number_neighbors);
        }

        void GEOM_GeometryMultigridQueue::IncrPriorityCV(unsigned long val_incr_point) 
        {
            /*--- Find the priority list ---*/
            short Number_Neighbors = Priority[val_incr_point];
            if (Number_Neighbors == -1) 
            {
                std::cout << "The CV " << val_incr_point << " is not in the priority list. (IncrPriorityCV)" << std::endl;
                exit(EXIT_FAILURE);
            }

            /*--- Remove the control volume ---*/
            RemoveCV(val_incr_point);

            /*--- Increase the priority ---*/
            AddCV(val_incr_point, Number_Neighbors + 1);
        }

        void GEOM_GeometryMultigridQueue::RedPriorityCV(unsigned long val_red_point)
        {
            /*--- Find the priority list ---*/
            short Number_Neighbors = Priority[val_red_point];
            if (Number_Neighbors == -1) 
            {
                std::cout << "The CV " << val_red_point << " is not in the priority list. (RedPriorityCV)" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (Number_Neighbors != 0) 
            {
                /*--- Remove the control volume ---*/
                RemoveCV(val_red_point);

                /*--- Increase the priority ---*/
                AddCV(val_red_point, Number_Neighbors - 1);
            }
        }

        void GEOM_GeometryMultigridQueue::VisualizeQueue(void) 
        {
            unsigned short iPoint;
            unsigned long jPoint;

            std::cout << std::endl;
            for (iPoint = 0; iPoint < QueueCV.size(); iPoint++) 
            {
                std::cout << "Number of neighbors " << iPoint << ": ";
                for (jPoint = 0; jPoint < QueueCV[iPoint].size(); jPoint++) 
                {
                    std::cout << QueueCV[iPoint][jPoint] << " ";
                }
                std::cout << std::endl;
            }
        }

        void GEOM_GeometryMultigridQueue::VisualizePriority(void)
        {
            unsigned long iPoint;

            for (iPoint = 0; iPoint < nPoint; iPoint++)
                std::cout << "Control Volume: " << iPoint << " Priority: " << Priority[iPoint] << std::endl;
        }

        long GEOM_GeometryMultigridQueue::NextCV(void)
        {
            if (QueueCV.size() != 0) return QueueCV[QueueCV.size() - 1][0];
            else return -1;
        }

        bool GEOM_GeometryMultigridQueue::EmptyQueue(void) 
        {
            unsigned short iPoint;

            /*--- In case there is only the no agglomerated elements,
            check if they can be agglomerated or we have already finished ---*/
            bool check = true;

            if (QueueCV.size() == 1) 
            {
                for (iPoint = 0; iPoint < QueueCV[0].size(); iPoint++) 
                {
                    if (RightCV[QueueCV[0][iPoint]]) { check = false; break; }
                }
            }
            else 
            {
                for (iPoint = 1; iPoint < QueueCV.size(); iPoint++)
                    if (QueueCV[iPoint].size() != 0) { check = false; break; }
            }

            return check;
        }

        unsigned long GEOM_GeometryMultigridQueue::TotalCV(void) 
        {
            unsigned short iPoint;
            unsigned long TotalCV;

            TotalCV = 0;
            for (iPoint = 0; iPoint < QueueCV.size(); iPoint++)
                if (QueueCV[iPoint].size() != 0) 
                { 
                    TotalCV += QueueCV[iPoint].size();
                }

            return TotalCV;
        }

        void GEOM_GeometryMultigridQueue::Update(unsigned long iPoint, GEOM_Geometry *fine_grid) 
        {
            unsigned short iNode;
            unsigned long jPoint;

            RemoveCV(iPoint);
            for (iNode = 0; iNode < fine_grid->node[iPoint]->GetnPoint(); iNode++)
            {
                jPoint = fine_grid->node[iPoint]->GetPoint(iNode);
                if (fine_grid->node[jPoint]->GetAgglomerate() == false)
                    IncrPriorityCV(jPoint);
            }
        }
    }
}