#ifndef __WINDOW_H
#define __WINDOW_H

#include <qgraphicsview.h>
#include <qpoint.h>
#include <qvector.h>

/* Let this class have access to objective values of the individuals
 * in the parent population.
 *
 * Once it has access to the double array then we create a QPointF
 * object with the x, y values (corresponding to objective 1 and 2)
 *
 * We then add that QPointF object to a vector and send it off to
 * the QGraphicsScene to show on the canvas*/

class Window: public QGraphicsView
{
  public:
    QGraphicsScene * scene;
    Window(QWidget *parent = 0);
    QVector <QPointF> points;
    
    void drawRectangle();
    
 };

#endif
