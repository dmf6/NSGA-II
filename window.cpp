#include "window.h"

#include <QtGui>

Window::Window(QWidget *parent) : QGraphicsView(parent) {
    scene = new QGraphicsScene();
    this->setScene(scene);
}

void Window::drawRectangle() {
    QPen blackPen(Qt::black);
    QBrush blueBrush(Qt::blue);
    blackPen.setWidth(6);
    QGraphicsEllipseItem *ellipse = scene->addEllipse(10,10, 200, 200, blackPen, blueBrush);
    ellipse->setFlag(QGraphicsItem::ItemIsMovable);
    
    // QGraphicsRectItem* item1 = new QGraphicsRectItem(0,0,100,100);
    // item1->setBrush(QBrush(Qt::red));
    // scene->addItem(item1);
}


    


