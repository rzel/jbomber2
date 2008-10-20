/*
 *  IMPORTANT: THIS FILE IS AUTOGENERATED
 *
 *  Any modifications to this file will be lost when regenerating it.
 *  Modify the corresponding .jpp file instead, and regenerate this file.
 */




import javax.swing.*;
import javax.swing.event.*;
import java.awt.event.*;
import java.awt.*;
import javax.swing.colorchooser.AbstractColorChooserPanel;
import java.util.Vector;




public class ColorSelector extends JDialog implements MouseListener, ChangeListener, ActionListener {
	public static final long serialVersionUID = 1;
	private Vector listeners = new Vector();

	public void addChangeListener(ChangeListener l) {
		listeners.add(l);
	}

	JColorChooser colorchooser;

	public ColorSelector(JFrame owner) {
		super(owner, "Color Chooser", true);
		setLayout(new FlowLayout());
		setResizable(false);
		addMouseListener(this);

		colorchooser = new JColorChooser();
		colorchooser.getSelectionModel().addChangeListener(this);


		AbstractColorChooserPanel[] panels = colorchooser.getChooserPanels();
		AbstractColorChooserPanel[] panels2 = new AbstractColorChooserPanel[panels.length + 1];
		panels2[0] = new ColorSelectorComponent();

		for (int i = 0; i < panels.length; ++i) {
			colorchooser.removeChooserPanel(panels[i]);
			panels2[i + 1] = panels[i];
		}

		colorchooser.setChooserPanels(panels2);


		add(colorchooser);

		pack();
	}

	public void stateChanged(ChangeEvent e) {
		if (e.getSource().equals(colorchooser.getSelectionModel())) {
			Color c = colorchooser.getColor();
			for (int i = 0 ; i < listeners.size(); ++i) {
				ChangeListener l = (ChangeListener)listeners.get(i);
				l.stateChanged(new ChangeEvent(c));
			}
		}
	}

	public void mousePressed(MouseEvent e) {
	}

	public void setVisible(boolean b) {
		super.setVisible(b);
	}

	void setSelectedColor(Color c) {
		colorchooser.getSelectionModel().setSelectedColor(c);
	}

	public void actionPerformed(ActionEvent e) {}
	public void mouseReleased(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseClicked(MouseEvent e) {}
}
