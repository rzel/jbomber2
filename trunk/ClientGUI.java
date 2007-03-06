import javax.swing.*;

import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.text.*;
import java.net.*;

/**
 * A <code>ServerGUI</code> object is a graphical interface to
 * start and stop a server, it also contains ways to keep track
 * of error messages returned by the server.
 */
public class ClientGUI extends JFrame
{
    // Something
    public static final long serialVersionUID = 1;

    // JTextArea to write server messages to
    private JTextArea messages = new JTextArea();

    // Panel to hold start_button and stop_button
    private JPanel top_panel = new JPanel();
    private JPanel bottom_panel = new JPanel();

    // Buttons to start, stop and save
    private JButton connect_button = new JButton("Connect");
    private JButton disconnect_button = new JButton("Disconnect");
    private JButton save_button = new JButton("Save log");

    // The NetworkClient that is used for this GUI
    private DelegatorClient client;

    // The LogPrintStream that is used by server to print messages
    private LogPrintStream log = new LogPrintStream(messages);

    /**
     * Creates a new <code>ClientGUI</code> object and initializes
     * all the GUI widgets.
     */
    public ClientGUI()
    {
        super("Bomberman - Client");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());

        // Initialize connect and disconnect buttons
        connect_button.addActionListener(new ConnectButtonHandler());
        disconnect_button.addActionListener(new DisconnectButtonHandler());
        disconnect_button.setEnabled(false);

        // Initialize save button
        save_button.addActionListener(new SaveButtonHandler());

        // Initialize top_panel
        top_panel.setLayout(new BorderLayout());
        top_panel.add(connect_button, BorderLayout.WEST);
        top_panel.add(disconnect_button, BorderLayout.EAST);

        // Initialize bottom_panel
        bottom_panel.setLayout(new FlowLayout(FlowLayout.CENTER));
        bottom_panel.add(save_button);

        messages.setEditable(false);

        // Add widgets to JFrame
        add(new JScrollPane(messages), BorderLayout.CENTER);
        add(top_panel, BorderLayout.NORTH);
        add(bottom_panel, BorderLayout.SOUTH);

        setSize(700, 900);

        setLocationRelativeTo(null);
        setVisible(true);
    }

    /**
     * A <code>StartButtonHandler</code> handles clicks on
     * a start button.
     */
    class ConnectButtonHandler implements ActionListener
    {
        /**
         * Handles click on a start button.
         */
        public void actionPerformed(ActionEvent e)
        {
            // Set button states
            connect_button.setEnabled(false);
            disconnect_button.setEnabled(true);

            // Try to start server
            try
            {
                InetAddress addr = InetAddress.getByName(Config.SERVER_ADDR);
                client = new DelegatorClient(addr,
                                             Config.SERVER_PORT,
                                             log,
                                             "Test");
            }
            catch (IOException ioe)
            {
                // If the server could not be started, give an error message.
                messages.append(ioe.getMessage() + "\n");
            }
        }
    }

    /**
     * A <code>StopButtonHandler</code> handles clicks on
     * a stop button.
     */
    class DisconnectButtonHandler implements ActionListener
    {
        /**
         * Handles click on a stop button.
         */
        public void actionPerformed(ActionEvent e)
        {
            // Stop server.
            client.close();
            client = null;
            // Set button states
            connect_button.setEnabled(true);
            disconnect_button.setEnabled(false);
        }
    }

    /**
     * A <code>SaveButtonHandler</code> handles clicks on
     * a save button.
     */
    class SaveButtonHandler implements ActionListener
    {
        /**
         * Handles click on a stop button.
         */
        public void actionPerformed(ActionEvent e)
        {
            // Save log contents to a file with name:
            // "Server_dd-mm-yyyy_hh.mm.ss.txt"

            // To format dates and times nicely with preceding 0's
            NumberFormat nf = NumberFormat.getInstance();
            nf.setMinimumIntegerDigits(2);

            // To get dates and times
            Calendar calendar = Calendar.getInstance();

            // Get date
            String day   = nf.format(calendar.get(calendar.DAY_OF_MONTH));
            String month = nf.format(calendar.get(calendar.MONTH) + 1);
            String year  = "" + calendar.get(calendar.YEAR);
            String date  = day + "-" + month + "-" + year;

            // Get time
            String hour  = nf.format(calendar.get(calendar.HOUR_OF_DAY));
            String min   = nf.format(calendar.get(calendar.MINUTE));
            String sec   = nf.format(calendar.get(calendar.SECOND));
            String time  = hour + "." + min + "." + sec;

            String log_name = "Client_" + date + "_" + time + ".txt";

            // Save contents
            try
            {
                FileWriter writer = new FileWriter(Config.LOG_PATH + log_name);
                messages.write(writer);
                writer.close();
                messages.append("LOG saved to: " +
                                Config.LOG_PATH + log_name +
                                "\n");
            }
            catch (IOException ioe)
            {
                // Error while saving, write error to log
                messages.append("Error while saving, message: " +
                                ioe.getMessage() + "\n");
            }

        }
    }

    /**
     * Start the ServerGUI.
     *
     * @param args Not used.
     */
    public static void main(String[] args)
    {
        new ClientGUI();
    }
}