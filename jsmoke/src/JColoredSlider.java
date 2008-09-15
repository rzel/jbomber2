/*
 *  IMPORTANT: THIS FILE IS AUTOGENERATED
 *
 *  Any modifications to this file will be lost when regenerating it.
 *  Modify the corresponding .jpp file instead, and regenerate this file.
 */




import java.awt.*;
import javax.swing.*;
import javax.swing.plaf.basic.*;
import java.util.*;
import java.awt.GradientPaint;

public class JColoredSlider extends JSlider
{
  private Color minColor;
  private Color maxColor;

  public JColoredSlider()
  {
    super();
  }

  public JColoredSlider(BoundedRangeModel brm)
  {
    super(brm);
  }

  public JColoredSlider(int orientation)
  {
    super(orientation);
  }

  public JColoredSlider(int min, int max)
  {
    super(min, max);
  }

  public JColoredSlider(int min, int max, int value)
  {
    super(min, max, value);
  }

  public JColoredSlider(int orientation, int min, int max, int value)
  {
    super(orientation, min, max, value);
  }

  public void setMinColor(Color color)
  {
    minColor = color;
  }

  public void setMaxColor(Color color)
  {
    maxColor = color;
  }

  public void paintComponent(Graphics gr)
  {
    float width = (float)getSize().getWidth();
    float height = (float)getSize().getHeight();

    Graphics2D g = (Graphics2D)gr;
    g.setColor(getBackground());
    g.fillRect(0,0,(int)width,(int)height);
    g.setPaint(new GradientPaint(0, 0, minColor, width, height, maxColor));
    g.fillRect(getX(),(int)height / 2,(int)width, (int)height);

    ((BasicSliderUI)getUI()).paintLabels(g);
    ((BasicSliderUI)getUI()).paintThumb(g);
  }
}

class JHueSlider extends JColoredSlider
{
  protected float saturation = 0f;
  protected float value = 100f;

  public JHueSlider()
  {
    super();
  }

  public JHueSlider(BoundedRangeModel brm)
  {
    super(brm);
  }

  public JHueSlider(int orientation)
  {
    super(orientation);
  }

  public JHueSlider(int min, int max)
  {
    super(min, max);
  }

  public JHueSlider(int min, int max, int value)
  {
    super(min, max, value);
  }

  public JHueSlider(int orientation, int min, int max, int value)
  {
    super(orientation, min, max, value);
  }

  public void setSaturationAndValue(float saturation, float value)
  {
    this.saturation = saturation;
    this.value = value;
  }

  public void paintComponent(Graphics gr)
  {
    float width = getWidth();
    float height = getHeight();

    Graphics2D g = (Graphics2D)gr;
    g.setColor(getBackground());
    g.fillRect(0,0,(int)width,(int)height);
    int precision = (int)width;

    for(int i = 0; i < precision; ++i)
    {
      float hue = i * ((float)ColorSelectorComponent.MAX_HUE /(float)precision);
      g.setColor(ColorConverter.hsvToColor(hue, saturation, value));
      g.fillRect(getX() + (int)((width / precision) * i),
                 (int)height / 2,
                 (int)((width / precision) * (i + 1)),
                 (int)height);
    }

    ((BasicSliderUI)getUI()).paintLabels(g);
    ((BasicSliderUI)getUI()).paintThumb(g);
  }
}

class JSliderlessSlider extends JColoredSlider
{
    private float[][] colors;
    private int banding;
    private int count;

  public JSliderlessSlider(BoundedRangeModel brm, float[][] colors, int count)
  {
    super();
    this.colors = colors;
    this.count  = count;
  }

  public JSliderlessSlider(BoundedRangeModel brm)
  {
    super(brm);
  }

  public JSliderlessSlider(int orientation)
  {
    super(orientation);
  }

  public JSliderlessSlider(int min, int max)
  {
    super(min, max);
  }

  public JSliderlessSlider(int min, int max, int value)
  {
    super(min, max, value);
  }

  public JSliderlessSlider(int orientation, int min, int max, int value)
  {
    super(orientation, min, max, value);
  }

  public void setColors(float[][] colors, int count) {
    this.colors = colors;
    repaint();
  }

  public void setCount(int count) {
      this.count = count;
      repaint();
  }

  public void paintComponent(Graphics gr)
  {
    double width = getWidth();
    double height = getHeight();

    Graphics2D g = (Graphics2D)gr;
    g.setColor(getBackground());
    g.fillRect(0,0,(int)width,(int)height);

    float k = 0;
    for(int i = 0; i < width; ++i) {
    	double pos = ((double)i) / width;
    	pos *= count ; pos = (int)(pos + .5); pos /= count;
    	int c = (int)(16384 * pos);
    	g.setColor(new Color(colors[c][0], colors[c][1], colors[c][2]));
    	g.drawLine(i, 0, i, (int)height);
    }




/*    for(int i = 0; i < count - 1; ++i) {
        Color color1 = new Color(colors[i][0], colors[i][1], colors[i][2]);
        Color color2 = new Color(colors[i + 1][0], colors[i + 1][1], colors[i + 1][2]);
        GradientPaint gp = new GradientPaint(k, 0, color1, k + width/count, height, color2);
        g.setPaint(gp);
        g.fillRect((int)k, 0, (int)k + getWidth()/count + 1, (int)height);
        k += getWidth()/count;
    }
*/

    g.setColor(Color.BLACK);
    g.drawRect(0, 0, (int)width - 1, (int)height - 1);
    ((BasicSliderUI)getUI()).paintLabels(g);
    ((BasicSliderUI)getUI()).paintTicks(g);
  }
}